import ee
ee.Authenticate()
ee.Initialize(project='')

import pandas as pd
import geopandas as gpd
import exactextract as eet
import pickle

from osgeo import gdal
from shapely import Polygon



##  Atlantic Canada

uplon = -52.05
uplat = 51.05
lowlon = -69.05
lowlat = 43.05

bound = [lowlon, lowlat, uplon, uplat]
atl_roi = ee.Geometry.Rectangle(bound)


##  Quebec City - Windsor corridor

uplon = -69.05
uplat = 48.20
lowlon = -86.05
lowlat = 41.20

bound = [lowlon, lowlat, uplon, uplat]
qw_roi = ee.Geometry.Rectangle(bound)



## Downloading from GEE

def gee_modis_lc_type3(roi, file_name, folder_name):

    landcover_collection = ee.ImageCollection('MODIS/061/MCD12Q1') \
                        .filterDate("2023-01-01", "2024-01-01") \
                        .map(lambda x: x.clip(roi)) \
                        .select(["LC_Type3"])

    ##  annually aggregated images
    landcover = landcover_collection.toList(landcover_collection.size())
    landcover = ee.Image(landcover.get(0))

    ##  coordinate matrix
    coordinates = landcover.pixelLonLat().select(['longitude','latitude'])

    task_1 = ee.batch.Export.image.toDrive(**{
            'image': landcover,
            'description': file_name,
            'folder': folder_name,
            'region': roi,
            'scale': 500,
            'fileFormat':'GeoTIFF',
            'crs':'EPSG:4326'
    })
    task_1.start()

    task_2 = ee.batch.Export.image.toDrive(**{
            'image': coordinates,
            'description': file_name + "_coordinates",
            'folder': folder_name,
            'region': roi,
            'scale': 500,
            'fileFormat':'GeoTIFF',
            'crs':'EPSG:4326'
    })
    task_2.start()



## Preload GEE exported TIF images

def preload_lc(file_name, folder_name):
    lc_rast = gdal.Open(folder_name + file_name)
    return lc_rast



## Preload TEMPO coordinate matrix (refer to "/download/tempo.py")

def preload_tempo_coords(file_name, folder_name):

    # arrays
    with open(folder_name + file_name, "rb") as f:
        coord = pickle.load(f)
    # geodataframe
    coord_gdf = gpd.GeoDataFrame({
        "longitude": coord[:,:,0].flatten(),
        "latitude": coord[:,:,1].flatten(),
    })
    # polygon geometry of 0.02 x 0.02 grids
    coord_gdf['geometry'] = coord_gdf.apply(lambda x: Polygon([ (x['longitude'] - 0.01, x['latitude'] + 0.01),
                                                                (x['longitude'] + 0.01, x['latitude'] + 0.01),
                                                                (x['longitude'] + 0.01, x['latitude'] - 0.01),
                                                                (x['longitude'] - 0.01, x['latitude'] - 0.01)
                                                            ]), axis=1)
    coord_gdf = coord_gdf.set_crs("epsg:4326")

    return coord_gdf



##  Resample 500 m x 500 m MODIS resolution to TEMPO grid resolution

def resample_lc(raster_lc, coord):

    # proportion of each land cover class in each TEMPO grid
    lc_prop = eet.exact_extract(raster_lc, coord, "frac")
    # presence of land classes in each TEMPO grid
    lc_unique_values = eet.exact_extract(raster_lc, coord, 'unique')

    # write to a new dataframe
    lc_resamp = pd.concat(
        [
            pd.DataFrame({
                u: v
                for u, v in zip(lc_unique_values[n]['properties']['unique'], lc_prop[n]['properties']['frac'])
            },
            index = [0])
            for n in range(coord.shape[0])
        ]
    )
    lc_resamp = lc_resamp.reset_index(drop=True)
    lc_resamp = lc_resamp.fillna(0)

    # rename columns with MODIS Land Cover Type 3 classes
    classes = [
        "water",
        "grassland",
        "shrubland",
        "cropland",
        "savannas",
        "forest_1",
        "forest_2",
        "forest_3",
        "forest_4",
        "barren",
        "urban"
    ]
    lc_resamp.columns = [classes[int(x)] for x in lc_resamp.columns]

    # combine classes
    lc_resamp['forest'] = lc_resamp[["forest_1","forest_2","forest_3","forest_4"]].sum(axis=1)
    lc_resamp['herbaceous'] = lc_resamp[["grassland","shrubland","cropland"]].sum(axis=1)

    lc_resamp = lc_resamp.drop(
        ["forest_1","forest_2","forest_3","forest_4",
         "grassland","shrubland","cropland"]
    )

    return lc_resamp



##  Reclassification rules

def land_cover_classifier(x):

    x = pd.DataFrame(x).transpose()
    max_cover = x.columns[np.argmax(x)]

    if np.max(x) >= 0.90:
        reclass = max_cover

    elif x['water'].values >= 0.30:
        if x['water'].values <= 0.60 and len(x[x >= 0.20]) >= 3 and max_cover != "water":
            reclass = "mixed with " + max_cover + "-water boundary"
        elif x['water'].values <= 0.60 and len(x[x >= 0.20]) >= 3 and max_cover == "water":
            reclass = "mixed with mainly water"
        elif x['water'].values <= 0.60 and len(x[x >= 0.20]) < 3:
            if x['urban'].values >= 0.20:
                reclass = "urban-water boundary"
            elif x['barren'].values >= 0.20:
                reclass = "barren-water boundary"
            elif x["forest"].values >= 0.20 or x["savannas"].values >= 0.20 or x["herbaceous"].values >= 0.20:
                reclass = "vegetated-water boundary"
            else:
                reclass = "mixed with mainly water"
        elif x['water'].values > 0.6 and len(x[x >= 0.20]) >= 2:
            if x['urban'].values >= 0.20:
                reclass = "urban-water boundary"
            elif x['barren'].values >= 0.20:
                reclass = "barren-water boundary"
            elif x["forest"].values >= 0.20 or x["savannas"].values >= 0.20 or x["herbaceous"].values >= 0.20:
                reclass = "vegetated-water boundary"
        else:
            reclass = "water"

    elif x['urban'].values >= 0.30:
        if x['urban'].values <= 0.60 and len(x[x >= 0.20]) >= 3 and max_cover != "urban":
            reclass = "mixed with " + max_cover + "-urban boundary"
        elif x['urban'].values <= 0.60 and len(x[x >= 0.20]) >= 3 and max_cover == "urban":
            reclass = "mixed with mainly urban"
        elif x['urban'].values <= 0.60 and len(x[x >= 0.20]) < 3:
            if x['water'].values >= 0.20:
                reclass = "urban-water boundary"
            elif x["barren"].values >= 0.20:
                reclass = "urban-barren boundary"
            elif x["forest"].values >= 0.20 or x["savannas"].values >= 0.20 or x["herbaceous"].values >= 0.20:
                reclass = "urban-vegetated boundary"
            else:
                reclass = "mixed with mainly urban"
        elif x['urban'].values > 0.6 and len(x[x >= 0.20]) >= 2:
            if x['water'].values >= 0.20:
                reclass = "urban-water boundary"
            elif x["barren"].values >= 0.20:
                reclass = "urban-barren boundary"
            elif x["forest"].values >= 0.20 or x["savannas"].values >= 0.20 or x["herbaceous"].values >= 0.20:
                reclass = "urban-vegetated boundary"
        else:
            reclass = "urban"

    elif x['barren'].values >= 0.30:
        if x['barren'].values <= 0.60 and len(x[x >= 0.20]) >= 3:
            if x['water'].values >= 0.20 or x['urban'].values >= 0.20:
                reclass = "mixed"
            else:
                reclass = "vegetated-barren boundary"
        elif x['barren'].values <= 0.6 and len(x[x >= 0.20]) < 3:
            if x['water'].values >= 0.20:
                reclass = 'barren-water boundary'
            elif x['urban'].values >= 0.20:
                reclass = 'urban-barren boundary'
            elif x["forest"].values >= 0.20 or x["savannas"].values >= 0.20 or x["herbaceous"].values >= 0.20:
                reclass = "vegetated-barren boundary"
            else:
                reclass = "mixed with mainly barren"
        elif x['barren'].values > 0.6 and len(x[x >= 0.20]) >= 2:
            if x['water'].values >= 0.20:
                reclass = 'barren-water boundary'
            elif x['urban'].values >= 0.20:
                reclass = 'urban-barren boundary'
            elif x["forest"].values >= 0.20 or x["savannas"].values >= 0.20 or x["herbaceous"].values >= 0.20:
                reclass = "vegetated-barren boundary"
        else:
            reclass = "barren"

    else:
        if len(x[x >= 0.30]) >= 2:
            reclass = "mixed vegetation"
        else:
            reclass = "mixed"

    return reclass

