import xarray as xr
import geopandas as gpd
import pandas as pd
import numpy as np
import pickle
import os
import time

import warnings
warnings.filterwarnings("ignore")

import earthaccess
earthaccess.login()


##  Atlantic Canada

boundary_lat_min = 43.05
boundary_lat_max = 51.05
boundary_lon_min = -69.05
boundary_lon_max = -52.05


##  Quebec City - Windsor corridor

boundary_lat_min = 41.20
boundary_lat_max = 48.20
boundary_lon_min = -86.05
boundary_lon_max = -69.05



##  spatial clipping function for tempo

def generate_geographical_subset(xarray, latmin, latmax, lonmin, lonmax):
    return xarray.where(
            (xarray.latitude < latmax) &
            (xarray.latitude > latmin) &
            (xarray.longitude < lonmax) &
            (xarray.longitude > lonmin),
            drop=True)



##  @ user-defined input paramters
##  @ 2024-09-25; only beta validation (V03) for Level 3 at that time
##  @ 2024-10-08; filter out unvalidated level 3 data (V01) with granule name string search 

def get_query(pollutant, level, lonmin, latmin, lonmax, latmax, time_start, time_end):
    
    product_name = "TEMPO_" + pollutant + "_" + level
    
    Query =  earthaccess.granule_query() \
                .short_name(product_name) \
                .readable_granule_name(["*V03*"]) \
                .cloud_hosted(True) \
                .bounding_box(lonmin, latmin, lonmax, latmax) \
                .temporal(time_start, time_end)
    
    cloud_granules = Query.get(Query.hits())

    return cloud_granules



##  monthly batches

year_st = [
    "2023-09-01",
    "2023-10-01",
    "2023-11-01",
    "2023-12-01",
    "2024-01-01",
    "2024-02-01",
    "2024-03-01",
    "2024-04-01",
    "2024-05-01",
    "2024-06-01",
    "2024-07-01",
    "2024-08-01",
]
year_ed = [
    "2023-09-30",
    "2023-10-31",
    "2023-11-30",
    "2023-12-31",
    "2024-01-31",
    "2024-02-29",
    "2024-03-31",
    "2024-04-30",
    "2024-05-31",
    "2024-06-30",
    "2024-07-31",
    "2024-08-31",
]

##  monthly batch downloading

def downloading(boundary_lat_min,
                boundary_lat_max,
                boundary_lon_min,
                boundary_lon_max):

    for s, e in zip(year_st, year_ed):
        
        ##  place query to earthdata.nasa.gov server
        cloud_granules = get_query(
            pollutant = "NO2",
            level = "L3",
            lonmin = boundary_lon_min,
            latmin = boundary_lat_min,
            lonmax = boundary_lon_max,
            latmax = boundary_lat_max,
            time_start = s,
            time_end = e
        )

        ##  initialize
        no2_collections = []
        no2_collections_support = []
        
        for n in range(len(cloud_granules)):

            files = earthaccess.download(cloud_granules[n], local_path="./Downloads/")
        
            ##  load the Latitude & Longitude coordinates
            ds_coord = xr.open_dataset(files[0])
        
            ##  load the data product features
            ds_product = xr.open_dataset(files[0], group='product')
            ds_support = xr.open_dataset(files[0], group='support_data')
        
            ##  assign the data product features with Latitude & Longitude coordinates
            ds_product = ds_product.assign_coords(ds_coord.coords)
            ds_support = ds_support.assign_coords(ds_coord.coords)
        
            ##  clip each dataset
            ds_coord = generate_geographical_subset(
                ds_coord,
                latmin = boundary_lat_min,
                latmax = boundary_lat_max,
                lonmin = boundary_lon_min,
                lonmax = boundary_lon_max
            )
            ds_product = generate_geographical_subset(
                ds_product,
                latmin = boundary_lat_min,
                latmax = boundary_lat_max,
                lonmin = boundary_lon_min,
                lonmax = boundary_lon_max
            )
            ds_support = generate_geographical_subset(
                ds_support,
                latmin = boundary_lat_min,
                latmax = boundary_lat_max,
                lonmin = boundary_lon_min,
                lonmax = boundary_lon_max
            )
            
            ## due to increasing order of latitude in .nc files, data arrays are inverted vertically
            ## axis 0 is time (dim=1), flipping at axis = 0 does not do the required operation, should be flipping axis = 1

            no2_collections.append({
                "Time": ds_coord.time.values,
                "Tropo_NO2": np.flip(ds_product['vertical_column_troposphere'].to_numpy(), axis=1), 
                "main_QA_flag": np.flip(ds_product['main_data_quality_flag'].to_numpy(), axis=1)
            })
            no2_collections_support.append({
                "Time": ds_coord.time.values,
                "eff_cloud_fraction": np.flip(ds_support['eff_cloud_fraction'].to_numpy(), axis=1),
                "snow_ice_fraction": np.flip(ds_support['snow_ice_fraction'].to_numpy(), axis=1),
                "amf_cloud_fraction": np.flip(ds_support['amf_cloud_fraction'].to_numpy(), axis=1),
            })
        
            ds_coord.close()
            ds_product.close()
            ds_support.close()
        
            os.remove("./Downloads/" + str(cloud_granules[n]).split("Data: ")[1].split("/")[-1].replace("']", ""))
        
            time.sleep(2.5)
        
            if (n + 1) % 25 == 0:
                print("completed: " + str(n+1))


        with open("./no2_" + s.replace("-01", "") + ".pkl", "wb") as f:
            pickle.dump(no2_collections_support, f)       
            
        time.sleep(5)



##  make a Lat/Lon coordinate matrix for the scene

def create_coord_matrix(xr_coord):

    coord_matrix = np.zeros((len(xr_coord.latitude), len(xr_coord.longitude), 2))

    for x in range(len(xr_coord.latitude)):
        coord_matrix[x, :, 0] = xr_coord['longitude'].values
    for y in range(len(xr_coord.longitude)):
        coord_matrix[:, y, 1] = xr_coord['latitude'].values[::-1]

    return coord_matrix
