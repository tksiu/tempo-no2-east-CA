import ee
ee.Authenticate()
ee.Initialize(project='')

import geopandas as gpd
import pandas as pd
import numpy as np
import seaborn
from matplotlib import pyplot as plt
from osgeo import gdal
from functools import reduce


def extract_agglomerate_boundary(
        agglomerate_list = [
            ["Halifax","Lake Echo","Still Water Lake","Brookside"],
            ["Cape Breton - Sydney","Sydney Mines","Glace Bay","New Waterford","Howie Centre","Eskasoni 3"],
            ["Saint John","Quispamsis - Rothesay","Wells"],
            ["Québec","Shannon","Stoneham","Sainte-Brigitte-de-Laval","Saint-Henri-de-Lévis","Château-Richer"],
        ]
):

    pc_shp = gpd.read_file(f"./Population Center/lpc_000b21a_e.shp")
    pc = pc_shp[pc_shp["PCNAME"].isin([y for x in agglomerate_list for y in x])]
    pc_geom_gjs = gpd.GeoDataFrame(pc["geometry"])
    pc_geom_gjs = pc_geom_gjs.set_crs(epsg=3347)
    pc_geom_gjs = pc_geom_gjs.to_crs(epsg=4326)
    pc_geom_gjs_bounds = pc_geom_gjs.bounds
    pc = pd.concat([pc, pc_geom_gjs_bounds], axis=1)

    for a in agglomerate_list:
        for a_len in range(1, len(a)):
            pc.loc[pc["PCNAME"] == a[a_len], "PCNAME"] = a[0]

    pc = pc.groupby(["PCNAME"]).agg({"minx": np.min,
                                    "miny": np.min,
                                    "maxx": np.max,
                                    "maxy": np.max}).reset_index()
    
    return pc


def gee_export_era5(pc):

    names = ["Cape Breton", "Halifax", "Quebec", "Saint John"]

    for p in range(pc.shape[0]):

        uplon = pc['maxx'].iloc[p]
        uplat = pc['maxy'].iloc[p]
        lowlon = pc['minx'].iloc[p]
        lowlat = pc['miny'].iloc[p]

        bound = [lowlon, lowlat, uplon, uplat]
        pc_roi = ee.Geometry.Rectangle(bound)

        for band in ['u_component_of_wind_10m', 'v_component_of_wind_10m']:

            era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR') \
                    .filterDate("2023-09-01", "2024-09-01") \
                    .select(band) \
                    .map(lambda x: x.clip(pc_roi)) \
                    .reduce(ee.Reducer.mean())

            tasks = ee.batch.Export.image.toDrive(**{
                    'image': era5,
                    'description': names[p] + "_" + band,
                    'folder': "era5_wind",
                    'region': pc_roi,
                    'scale': 2500,
                    'fileFormat':'GeoTIFF',
                    'crs':'EPSG:4326'
            })
            tasks.start()


def compile_information(city_names, folder):

    features = {}
    features_df = []

    for p in city_names:

        features[p] = []
        feature_coords = gdal.Open(folder + p + "_coords.tif")
        lat = np.array(feature_coords.GetRasterBand(2).ReadAsArray())
        lon = np.array(feature_coords.GetRasterBand(1).ReadAsArray())

        for band in ['u_component_of_wind_10m', 'v_component_of_wind_10m']:

            tem = gdal.Open(folder + p + "_" + band + ".tif")
            tem = np.array(tem.GetRasterBand(1).ReadAsArray())
            obs = pd.DataFrame({
                "longitude": lon.flatten(),
                "latitude": lat.flatten(),
                band: tem.flatten()
            })

            features[p].append(obs)

        tem = reduce(lambda df1, df2: pd.merge(df1, df2, on=["longitude", "latitude"]), features[p])
        tem = tem.dropna()
        tem['pname'] = p

        features_df.append(tem)

        tem = reduce(lambda df1, df2: pd.merge(df1, df2, on=["longitude", "latitude"]), features[p])
        tem = tem.dropna()
        tem['pname'] = p

        features_df.append(tem)

    features_df = pd.concat(features_df)

    return features_df


def plot_annual_wind_components(plot_df):

    plot_df['pname'] = plot_df['pname'].apply(
        lambda x: "Quebec City" if x == "Quebec" else x
    )

    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    seaborn.scatterplot(
        data = plot_df,
        x = 'u_component_of_wind_10m',
        y = 'v_component_of_wind_10m',
        hue = "pname", style = "pname",
        s = 120,
        palette = "Dark2",
        ax = ax
    )

    ax.set_ylim([np.max(plot_df['v_component_of_wind_10m']) * -1.1, np.max(plot_df['v_component_of_wind_10m']) * 1.1])
    ax.set_xlim([np.max(plot_df['u_component_of_wind_10m']) * -0.8, np.max(plot_df['u_component_of_wind_10m']) * 1.1])
    ax.axhline(0, color='black', linewidth=.5)
    ax.axvline(0, color='black', linewidth=.5)
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.set_xlabel("UGRD")
    ax.set_ylabel("VGRD")
    ax.legend(title="", fontsize=14)
