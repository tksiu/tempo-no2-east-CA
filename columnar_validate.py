import numpy as np
import pandas as pd
import geopandas as gpd
import rioxarray as rxr
import xarray as xr
import rasterio as rio

import os
import pickle
import datetime
import seaborn
import scipy

from matplotlib import pyplot as plt



def read_tempo_time():

    # TEMPO: download process referred to "/download/tempo.py" and "qa_aggregate.py"
    tempo_folder = "./TEMPO/"
    input_subfolders = [tempo_folder + "Source/SON/",
                        tempo_folder + "Source/DJF/",
                        tempo_folder + "Source/MAM/",
                        tempo_folder + "Source/JJA/"]

    tempo_month_time = []

    for n in range(len(input_subfolders)):
        
        no2_obs = [x for x in os.listdir(input_subfolders[n])
                     if "_support_data" not in x and "agg" not in x and "count" not in x]

        for obs in no2_obs: 

            with open(input_subfolders[n] + obs, "rb") as f:
                data = pickle.load(f)
                
            tempo_timestamp = [pd.to_datetime(x["Time"][0]) for x in data]
            tempo_month_time.append(tempo_timestamp)



def read_tropomi_time():

    # TROPOMI: download process referred to "/download/tropomi.py"
    tropomi_folder = "./TROPOMI/"
    tropomi_files = os.listdir(tropomi_folder)

    tropomi_timestamp = [x.replace("TROPOMI_NO2_", "").replace(".tif", "").replace("T", " ").replace("_", ":") 
                         for x in tropomi_files]

    # restrict +/- 15 minute difference from TEMPO
    tropomi_upbound = [datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") + datetime.timedelta(minutes=15)
                    for x in tropomi_timestamp]
    tropomi_lowbound = [datetime.datetime.strptime(x, "%Y-%m-%d %H:%M:%S") - datetime.timedelta(minutes=15)
                        for x in tropomi_timestamp]
    


def temporal_matching(tempo_month_time, tropomi_upbound, tropomi_lowbound):
    
    ##  retrieve all pairs (indices) of temporally matched (+/- 15 min) TEMPO and TROPOMI images
    match_idx = []

    for n in range(len(tempo_month_time)):

        match_idx_sub = []

        for m in range(len(tempo_month_time[n])):

            # restrict +/- 15 minute difference between TEMPO and TROPOMI
            pull_idx = [(m, t)
                        for t in range(len(tropomi_timestamp))
                        if tempo_month_time[n][m] < tropomi_upbound[t] and tempo_month_time[n][m] > tropomi_lowbound[t]]
            
            if len(pull_idx) > 0:
                match_idx_sub.append(pull_idx)
                
        match_idx.append(match_idx_sub)

    return match_idx



def read_tempo_coords(file_name="tempo_coord_array.pkl", folder_name="./"):

    ## Preload TEMPO coordinate matrix (refer to "/download/tempo.py")
    with open(folder_name + file_name, "rb") as f:
        coord = pickle.load(f)
        
    coord_gdf = gpd.GeoDataFrame({
        "longitude": coord[:,:,0].flatten(),
        "latitude": coord[:,:,1].flatten(),
    })

    tempo_coords_xr = xr.Dataset.from_dataframe(coord_gdf.set_index(['longitude','latitude']))
    tempo_coords_xr = tempo_coords_xr.rio.write_crs("epsg:4326")
    tempo_coords_xr = tempo_coords_xr.rename({"longitude": "x", "latitude": "y"})
    tempo_coords_xr = tempo_coords_xr.reindex(y=list(reversed(tempo_coords_xr.y)))



def read_tropomi_coords(file_name="tropomi_coord_array.pkl", folder_name="./"):

    ## Preload TROPOMI coordinate matrix (refer to "/download/tropomi.py")
    with open(folder_name + file_name, "rb") as f:
        tropomi_coords = pickle.load(f)

    tropomi_coords_gdf = pd.DataFrame({
            "longitude": tropomi_coords[:,:,0].flatten(),
            "latitude": tropomi_coords[:,:,1].flatten()
    })



def colocate_tempo_tropomi(data,                ## original image to be resampled:  TROPOMI
                           coords,              ## original image's grid coordinates
                           target_grid,         ## targeted grid coordinates for resampling:  TEMPO
                           resampling_scheme    ## resampling method
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



def operate_pull_2km(tempo_img,                                             ## TEMPO image
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
                     snow_fraction=0.01):                                   ## snow/ice fraction cutoff if applying QA screening
        
    """ 1) TEMPO images """

    tempo_tem = copy.deepcopy(tempo_img[main_col][0,:,:])

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
            "latitude": coord[:,:,1].flatten(),
            "longitude": coord[:,:,0].flatten(),
            "tempo_px": tempo_tem.flatten()
    })

    """ 2) TROPOMI images """

    tropomi_img = rxr.open_rasterio(tropomi_folder + tropomi_files[match_idx[month_idx][n][0][1]]).values
    tropomi_img = tropomi_img[0,:,:]

    # convert mol/m^2 to molec/cm^2
    tropomi_img = tropomi_img * 6.022e23 / 1e4

    tropomi_img = colocate_tempo_tropomi(tropomi_img, tropomi_coords_gdf, tempo_coords_xr, resampling_scheme)

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

    read_tempo_time()
    read_tropomi_time()
    read_tempo_coords()
    read_tropomi_coords()

    for folder in input_folders:

        no2_obs = [x for x in os.listdir(folder) if "_support_data" not in x and "no2" in x and "agg" not in x and "count" not in x] 
        no2_support = [x for x in os.listdir(folder) if "_support_data" in x and "no2" in x and "agg" not in x and "count" not in x]
        
        for obs, aux in zip(no2_obs, no2_support):    
            
            if len(match_idx[month_idx]) > 0:
                
                with open(folder + obs, "rb") as f:
                    tempo_collection = pickle.load(f)
        
                with open(folder + aux, "rb") as f:
                    tempo_collection_aux = pickle.load(f)
                    
                for n in range(len(match_idx[month_idx])):
                    
                    tempo_case = tempo_collection[match_idx[month_idx][n][0][0]]
                    tempo_aux = tempo_collection_aux[match_idx[month_idx][n][0][0]]
                    
                    df, arr = operate_pull_2km(
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
                    tropomi_match_time.append(tropomi_timestamp[match_idx[month_idx][n][0][1]])
                    
            month_idx += 1
            print("completed.")


    ###  Analysis (correlation)  ###
    
    # remove observations < 1e14 (extremely low values, only less than 1% of data)
    for x in range(len(tempo_data)):
        tempo_data[x][tempo_data[x] < 1e14] = np.nan
        tropomi_data[x][tropomi_data[x] < 1e14] = np.nan

    tempo_agg = np.nanmean(np.concatenate([x.reshape(1, x.shape[0], x.shape[1]) for x in tempo_data]), axis=0)
    tropomi_agg = np.nanmean(np.concatenate([x.reshape(1, x.shape[0], x.shape[1]) for x in tropomi_data]), axis=0)

    # add MODIS land cover classes (refer to "/download/modis_lc.py")
    tempo_lc = pd.read_csv("./lc.csv")

    # gather TEMPO and TROPOMI values for correlation analysis, stratified by land covers
    corr_df = pd.DataFrame({
        "tempo": tempo_agg.flatten(),
        "tropomi": tropomi_agg.flatten(),
        "lc": tempo_lc['reclass'].values
    })
    corr_df = corr_df.dropna()


    ###  Visualization  ###

    # scatter plot
    seaborn.scatterplot(data=corr_df, x="tempo", y="tropomi")

    # simple linear analysis
    sm.OLS(corr_df['tropomi'], sm.add_constant(corr_df['tempo']))\
        .fit()\
        .summary()
    
    # scatter plot with density infromation
    seaborn.regplot(data=corr_df, x="tempo", y="tropomi",
                    marker='o', line_kws={"color": "red"}, scatter_kws={'s':2, "color": "black"})
    seaborn.kdeplot(data=corr_df, x="tempo", y="tropomi", 
                    fill=True)
    plt.xlabel("TEMPO", fontsize=16)
    plt.ylabel("TROPOMI", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.text(np.nanmax(corr_df["tempo"])*1.05, np.nanmax(corr_df["tropomi"])*1.05, 
            "Spearman's R = {:.2f}".format(scipy.stats.spearmanr(corr_df.tempo, corr_df.tropomi)[0]),
            fontsize=16)

    # residual plot of linear regression
    seaborn.residplot(data=corr_df, x="tempo", y="tropomi", 
                      scatter_kws={'s':2, "color": "black"})
    plt.xlabel("TEMPO", fontsize=16)
    plt.ylabel("Residual", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    
    # pixelwise difference
    tempo_coords_xr['tempo_agg'] = (('y', 'x'), tempo_agg)
    tempo_coords_xr['tropomi_agg'] = (('y', 'x'), tropomi_agg)

    tempo_coords_xr = tempo_coords_xr.interpolate_na(dim="y", method="linear").interpolate_na(dim="x", method="linear")

    diff = tempo_coords_xr.tempo_agg - tempo_coords_xr.tropomi_agg
