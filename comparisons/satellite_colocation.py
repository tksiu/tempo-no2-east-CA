import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import xarray as xr
import rasterio as rio

import os
import copy
import pickle
import datetime




class temporal_matching:

    def __init__(self, tempo_folder: str, tropomi_folder: str):
        
        # TEMPO: download process referred to "/download/tempo.py"
        self.tempo_folder = tempo_folder
        self.tempo_subfolders = os.listdir(self.tempo_folder)
        
        # TROPOMI: download process referred to "/download/tropomi.py"
        self.tropomi_folder = tropomi_folder
        self.tropomi_files = os.listdir(self.tropomi_folder)


    def read_tempo_time(self):

        tempo_month_time = []

        for n in range(len(self.tempo_subfolders)):
            
            no2_obs = [x for x in os.listdir(self.tempo_subfolders[n])
                        if "_support_data" not in x and "agg" not in x and "count" not in x]

            for obs in no2_obs: 

                with open(self.tempo_subfolders[n] + obs, "rb") as f:
                    data = pickle.load(f)
                    
                tempo_timestamp = [pd.to_datetime(x["Time"][0]) for x in data]
                tempo_month_time.append(tempo_timestamp)
        
        return tempo_month_time

    def read_tropomi_time(self):

        tropomi_timestamp = [x.replace("TROPOMI_NO2_", "").replace(".tif", "").replace("T", " ").replace("_", ":") 
                            for x in self.tropomi_files]

        # restrict +/- 15 minute difference from TEMPO
        tropomi_upbound = [datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") + datetime.timedelta(minutes=15)
                        for x in tropomi_timestamp]
        tropomi_lowbound = [datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") - datetime.timedelta(minutes=15)
                            for x in tropomi_timestamp]
        
        return tropomi_timestamp, tropomi_upbound, tropomi_lowbound

    def get_match_index(self):

        self.tempo_month_time = self.read_tempo_time()
        self.tropomi_timestamp, self.tropomi_upbound, self.tropomi_lowbound = self.read_tropomi_time()
        
        ##  retrieve all pairs (indices) of temporally matched (+/- 15 min) TEMPO and TROPOMI images
        match_idx = []

        for n in range(len(self.tempo_month_time)):

            match_idx_sub = []

            for m in range(len(self.tempo_month_time[n])):

                # restrict +/- 15 minute difference between TEMPO and TROPOMI
                pull_idx = [(m, t)
                            for t in range(len(self.tropomi_timestamp))
                            if self.tempo_month_time[n][m] < self.tropomi_upbound[t] and self.tempo_month_time[n][m] > self.tropomi_lowbound[t]]
                
                if len(pull_idx) > 0:
                    match_idx_sub.append(pull_idx)
                    
            match_idx.append(match_idx_sub)

        return match_idx




class spatial_matching:

    def __init__(self, 
                 tempo_coord_folder: str, 
                 tempo_coord_file: str,
                 tropomi_coord_folder: str,
                 tropomi_coord_file: str):

        self.tempo_coords_xr, self.tempo_coords_gdf = self.read_coords(file_name = tempo_coord_file, folder_name = tempo_coord_folder)
        self.tropomi_coords_xr, self.tropomi_coords_gdf = self.read_coords(file_name = tropomi_coord_file, folder_name = tropomi_coord_folder)


    def read_coords(self, file_name="tempo_coord_array.pkl", folder_name="./"):

        ## Preload TEMPO coordinate matrix (refer to "/download/tempo.py")
        with open(folder_name + file_name, "rb") as f:
            coord = pickle.load(f)
            
        coord_gdf = gpd.GeoDataFrame({
            "longitude": coord[:,:,0].flatten(),
            "latitude": coord[:,:,1].flatten(),
        })

        coords_xr = xr.Dataset.from_dataframe(coord_gdf.set_index(['longitude','latitude']))
        coords_xr = coords_xr.rio.write_crs("epsg:4326")
        coords_xr = coords_xr.rename({"longitude": "x", "latitude": "y"})
        coords_xr = coords_xr.reindex(y=list(reversed(coords_xr.y)))

        return coords_xr, coord_gdf

    def colocate_tempo_tropomi(self,
                               data,                ## original image to be resampled:  TROPOMI
                               coords,                 ## original image's grid coordinates
                               target_grid,            ## targeted grid coordinates for resampling:  TEMPO
                               resampling_scheme       ## resampling method
                               ):

        df = pd.DataFrame({
            "latitude": coords.latitude,
            "longitude": coords.longitude,
            "observation": data.flatten()
        })

        df = xr.Dataset.from_dataframe(df.set_index(['longitude','latitude']))
        df = df.rio.write_crs("epsg:4326")
        df = df.rename({"longitude": "x", "latitude": "y"})
        df = df.transpose("y", "x")
        df = df.reindex(y=list(reversed(df.y)))

        df = df.rio.reproject_match(
                        target_grid,
                        resampling = resampling_scheme
        )
        return df

    def operate_2km(self,
                    tempo_img,                                             ## TEMPO image
                    tropomi_folder,                                        ## folder for TROPOMI series of images
                    match_idx,                                             ## temporally matched indices
                    month_idx,                                             ## iterate each month
                    main_col="Tropo_NO2",                                  ## TEMPO Tropospheric NO2 VCD field name
                    qa_flag_col=None,                                      ## specify QA flag field name if applying QA screening 
                    resampling_scheme=rio.enums.Resampling.bilinear,       ## resampling method
                    tempo_cloud=None,                                      ## specify TEMPO support data (cloud fraction) file list
                    cloud_col='amf_cloud_fraction',                        ## cloud fraction variable if applying QA screening
                    cloud_fraction=0.50,                                   ## cloud fraction cutoff if applying QA screening
                    tempo_snow=None,                                       ## specify TEMPO support data (snow fraction) file list
                    snow_col='snow_ice_fraction',                          ## snow/ice fraction variable if applying QA screening
                    snow_fraction=0.01                                     ## snow/ice fraction cutoff if applying QA screening
                    ):                                   
            
        """ 1) TEMPO images """

        tempo_tem = copy.deepcopy(tempo_img[main_col][0,:,:])
        tropomi_files = os.listdir(tropomi_folder)

        # screening with QA flag
        if qa_flag_col != None:
            screen = copy.deepcopy(tempo_img[qa_flag_col][0,:,:])
            tempo_tem[screen != 0] = np.nan

        # cloud cover < 0.50  (align with TROPOMI Level 3)
        if tempo_cloud != None:
            cloud_screen = copy.deepcopy(tempo_cloud[cloud_col][0,:,:])
            tempo_tem[cloud_screen > cloud_fraction] = np.nan
            
        # snow/ice cover < 0.01  (align with TROPOMI Level 3)
        if tempo_snow != None:
            snow_screen = copy.deepcopy(tempo_snow[snow_col][0,:,:])
            tempo_tem[snow_screen > snow_fraction] = np.nan

        # inverting along horizontal axis
        tempo_tem = np.flip(tempo_tem, axis=0)

        tempo_df = pd.DataFrame({
                "latitude": self.tempo_coords_xr[:,:,1].flatten(),
                "longitude": self.tempo_coords_xr[:,:,0].flatten(),
                "tempo_px": tempo_tem.flatten()
        })

        """ 2) TROPOMI images """

        tropomi_img = rxr.open_rasterio(tropomi_folder + tropomi_files[match_idx[month_idx][n][0][1]]).values
        tropomi_img = tropomi_img[0,:,:]

        # convert mol/m^2 to molec/cm^2
        tropomi_img = tropomi_img * 6.022e23 / 1e4

        tropomi_img = self.colocate_tempo_tropomi(tropomi_img, self.tropomi_coords_gdf, self.tempo_coords_xr, resampling_scheme)

        tropomi_df = tropomi_img.to_dataframe().reset_index()


        tempo_df = tempo_df[~pd.isnull(tempo_df["tempo_px"])]
        tempo_df = tempo_df.reset_index(drop=True)

        tropomi_df = tropomi_df[~pd.isnull(tropomi_df["observation"])]
        tropomi_df = tropomi_df.reset_index(drop=True)

        tropomi_df.columns = ["longitude","latitude","spatial_ref","tropomi_px"]

        df = tempo_df.merge(tropomi_df, on=["longitude","latitude"], how="inner")
            
        return df, (tempo_tem, tropomi_img.observation.values)





if __name__ == "main":

    ###  Extraction  ###

    data = []
    tempo_data = []
    tropomi_data = []

    tropomi_match_time = []
    tempo_match_time = []

    month_idx = 0
    tempo_folder = "./TEMPO/"               ## atlantic or QW-corridor specific
    tropomi_folder = "./TROPOMI/"           ## atlantic or QW-corridor specific

    input_folders = [                       ##  Example folder structure:  monthly collections of images (.pkl) by season
        tempo_folder + "Source/SON/",
        tempo_folder + "Source/DJF/",
        tempo_folder + "Source/MAM/",
        tempo_folder + "Source/JJA/"
    ]
    
    match_instance = temporal_matching(tempo_folder, tropomi_folder)
    match_idx = match_instance.get_match_index()

    colocation_instance = spatial_matching(tempo_coord_folder = "./coordinates/", 
                                            tempo_coord_file = "tempo_coord_array.pkl",         ## atlantic or QW-corridor specific
                                            tropomi_coord_folder = "./coordinates/",
                                            tropomi_coord_file = "tropomi_coord_array.pkl")     ## atlantic or QW-corridor specific

    for folder in input_folders:

        no2_obs = [x for x in os.listdir(folder) if "_support_data" not in x and "no2" in x and "agg" not in x and "count" not in x] 
        no2_support = [x for x in os.listdir(folder) if "_support_data" in x and "no2" in x and "agg" not in x and "count" not in x]
        
        for obs, aux in zip(no2_obs, no2_support):    
            
            if len(match_idx[month_idx]) > 0:
                
                ##  TEMPO main NO2 product
                with open(folder + obs, "rb") as f:
                    tempo_collection = pickle.load(f)

                ##  TEMPO supporting data for QA/QC
                with open(folder + aux, "rb") as f:
                    tempo_collection_aux = pickle.load(f)
                    
                for n in range(len(match_idx[month_idx])):
                    
                    tempo_case = tempo_collection[match_idx[month_idx][n][0][0]]
                    tempo_aux = tempo_collection_aux[match_idx[month_idx][n][0][0]]
                    
                    df, arr = colocation_instance.operate_2km(
                        tempo_case, 
                        match_idx, 
                        month_idx,
                        main_col="Tropo_NO2", 
                        qa_flag_col="main_QA_flag",
                        resampling_scheme = rio.enums.Resampling.average,
                        tempo_cloud=tempo_aux, 
                        cloud_col='amf_cloud_fraction', 
                        cloud_fraction=0.50,
                        tempo_snow=tempo_aux, 
                        snow_col='snow_ice_fraction', 
                        snow_fraction=0.01,
                    )
                    
                    data.append(df)
                    tempo_data.append(arr[0])
                    tropomi_data.append(arr[1])
                    
                    tempo_match_time.append(str(tempo_case["Time"])[2:-9].replace("T", " "))
                    tropomi_match_time.append(match_instance.tropomi_timestamp[match_idx[month_idx][n][0][1]])
                    
            month_idx += 1
            print("completed.")

