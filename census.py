import geopandas as gpd
import pandas as pd
import numpy as np
import shapely
import math


###  Data Source: 2021 Canadian Census of Population @ https://www12.statcan.gc.ca/census-recensement/2021/dp-pd/prof/index.cfm?Lang=E
###  Boundary Source:  Statistics Canada @ https://www12.statcan.gc.ca/census-recensement/2021/geo/sip-pis/boundary-limites/index2021-eng.cfm?year=21


# Initialization:  selected areas with high annual NO2 and large populations

agglomerate_list = [
    ["Halifax","Lake Echo","Still Water Lake","Brookside"],
    ["Cape Breton - Sydney","Sydney Mines","Glace Bay","New Waterford","Howie Centre","Eskasoni 3"],
    ["Saint John","Quispamsis - Rothesay","Wells"],
    ["Montréal","Châteauguay","Les Cèdres","Saint-Jérôme","Saint-Jean-sur-Richelieu","Saint-Sulpice","Beloei"],
    ["Québec","Shannon","Stoneham","Sainte-Brigitte-de-Laval","Saint-Henri-de-Lévis","Château-Richer"],
    ["Toronto","Oshawa","Milton","Georgetown","Mount Albert"],
    ["Windsor","Amherstburg"],
    ["Sarnia","Corunna","Petrolia","Wyoming"]
]



###  1) Population Centers for targeted cities

def pop_center():

    pc_shp_path = f"./Shapefiles/Population Center/"
    pc_shp = gpd.read_file(pc_shp_path + "lpc_000b21a_e.shp")

    # excluding irrelevant population centers

    pc = pc_shp[pc_shp["PCNAME"].isin([y for x in agglomerate_list for y in x])]
    pc = pc[~((pc["PCNAME"] == "Windsor") & (pc["LANDAREA"] < 20))]

    # extracting bounding boxes

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



def pc_interpolate_grids():

    ###  for intra-populaation-center analyses, downsample TEMPO 0.02 x 0.02 grids to 0.002 x 0.002 (200 m x 200 m) resolution

    pc_coords = {}
    pc_grids_coords = {}
    aspect = {}

    for n in range(pc.shape[0]):

        name = pc["PCNAME"].iloc[n]

        ##  calculate scene widths and heights (so as the aspect ratios) for selected PCs
        width = np.array((pc.maxx.iloc[n] - pc.minx.iloc[n]) / 0.002)
        height = np.array((pc.maxy.iloc[n] - pc.miny.iloc[n]) / 0.002)

        ##  consider buffering two grid-line edges for each dimension
        dim = (math.ceil(width), math.ceil(height))

        ##  starting from northend and westend, get grid corner and center coordinates 
        top_left_corner_h = np.array([pc.minx.iloc[n] + 0.002 * z for z in range(dim[0])]).reshape(-1).tolist()
        top_right_corner_h = top_left_corner_h[1:] + [top_left_corner_h[-1] + 0.002]
        top_left_corner_v = np.array([pc.maxy.iloc[n] - 0.002 * z for z in range(dim[1])]).reshape(-1).tolist()
        bottom_left_corner_v = top_left_corner_v[1:] + [top_left_corner_v[-1] - 0.002]
        centroids_h = top_left_corner_h + [top_right_corner_h[-1]]
        centroids_v = top_left_corner_v + [bottom_left_corner_v[-1]]
        centroids_h = [0.5 * (centroids_h[x] + centroids_h[x+1]) for x in range(len(centroids_h) - 1)]
        centroids_v = [0.5 * (centroids_v[x] + centroids_v[x+1]) for x in range(len(centroids_v) - 1)]

        ##  save the interpolated grids
        polys_coords = list()
        polys_ids = list()

        for i in range(dim[0]):
            for j in range(dim[1]):

                polys_coords.append(
                    shapely.Polygon(
                        [[top_left_corner_h[i], top_left_corner_v[j]],
                        [top_right_corner_h[i], top_left_corner_v[j]],
                        [top_right_corner_h[i], bottom_left_corner_v[j]],
                        [top_left_corner_h[i], bottom_left_corner_v[j]],
                        [top_left_corner_h[i], top_left_corner_v[j]],
                    ])
                )

                polys_ids.append(str(i) + "_" + str(j))

        grids_coords = gpd.GeoDataFrame({
            "id": polys_ids,
            "geometry": polys_coords
        })
        grids_coords['lon'] = grids_coords['geometry'].centroid.x
        grids_coords['lat'] = grids_coords['geometry'].centroid.y
        grids_coords['lon'] = grids_coords['lon'].apply(lambda x: round(x, 6))
        grids_coords['lat'] = grids_coords['lat'].apply(lambda x: round(x, 6))

        coords = [
            centroids_h,
            centroids_v
        ]

        pc_coords[name] = coords
        pc_grids_coords[name] = grids_coords
        aspect[name] = dim



###  2) Census Tracts for targeted cities (except Cape Breton)

def tracts(predicate):

    ct_shp_path = f"./Shapefiles/Census Tract/"
    ct_polys = gpd.read_file(ct_shp_path + "lct_000b21a_e.shp")

    ct_polys = ct_polys.set_crs(epsg=3347)
    ct_polys = ct_polys.to_crs(epsg=4326)

    # extracting all census tracts within or intersecting the outer bounding boxes of selected population centers

    assert predicte in ["within", "intersect"], "available options are: ['within', 'intersect']"

    collect_cts = {}

    for p in range(pc.shape[0]):

        pc_name = pc["PCNAME"].iloc[p]
        
        if predicte == "within":

            extract = ct_polys.within(
                shapely.Polygon([   
                    (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02),
                    (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                    (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                    (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02)
                ])
            )

        elif predicte == "intersect":

            extract = ct_polys.intersection(
                shapely.Polygon([   
                    (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02),
                    (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                    (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                    (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                    pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02)
                ])
            )
            
        extract = pd.DataFrame(extract).reset_index(drop=True)
        extract.columns = ["included"]

        extract["DGUID"] = ct_polys["DGUID"]
        extract["geometry"] = ct_polys["geometry"]
        extract["PCNAME"] = pc_name
        extract = extract[extract["included"] == True]

        extract = gpd.GeoDataFrame(extract).reset_index(drop=True)
        extract = extract.set_geometry("geometry")
        
        collect_cts[pc_name] = extract

    return collect_cts



### 3) Census Dissemination Areas (Cape Breton)

def dissemination_areas(predicate):

    da_shp_path = f"./Shapefiles/Dissemination Area/"
    da_polys = gpd.read_file(da_shp_path + "lda_000b21a_e.shp")

    da_polys = da_polys.set_crs(epsg=3347)
    da_polys = da_polys.to_crs(epsg=4326)

    # extracting all census dissemination areas within the outer bounding boxes of Cape Breton (no tracts available)

    assert predicte in ["within", "intersect"], "available options are: ['within', 'intersect']"

    collect_das = {}

    pc_name = "Cape Breton - Sydney"

    if predicte == "within":

        extract = da_polys.within(
            shapely.Polygon([   
                (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02),
                (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02)
            ])
        )

    elif predicte == "intersect":

        extract = da_polys.intersection(
            shapely.Polygon([   
                (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02),
                (pc[pc["PCNAME"] == pc_name].values[0][1] - 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][4] + 0.02),
                (pc[pc["PCNAME"] == pc_name].values[0][3] + 0.02,
                pc[pc["PCNAME"] == pc_name].values[0][2] - 0.02)
            ])
        )

    extract = pd.DataFrame(extract).reset_index(drop=True)
    extract.columns = ["included"]

    extract["DGUID"] = da_polys["DGUID"]
    extract["geometry"] = da_polys["geometry"]
    extract["PCNAME"] = pc_name
    extract = extract[extract["included"] == True]
    
    extract = gpd.GeoDataFrame(extract).reset_index(drop=True)
    extract = extract.set_geometry("geometry")
    
    collect_das[pc_name] = extract

    return collect_das
