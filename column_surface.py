import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr

import os
import pickle
import shapely
import datetime
import seaborn

from matplotlib import pyplot as plt
from numpy import dot
from numpy.linalg import norm
from itertools import chain



def load_monitor_obs():

    #  read AirNow observations (refer to "/download/airnow.py" for details; combine all files into one)
    valid_monitor_obs = pd.read_excel(f"/AirNow/NO2/airnow_no2.xlsx")

    #  re-sorting
    valid_monitor_obs = valid_monitor_obs.sort_values(['LONGITUDE(deg)', 'LATITUDE(deg)','Timestamp(UTC)']).reset_index(drop=True)

    # create point (vector) dataframe for the regulatory monitor locations
    valid_pts = valid_monitor_obs[['LONGITUDE(deg)', 'LATITUDE(deg)']].drop_duplicates()
    valid_pts['geometry'] = valid_pts.apply(lambda x: shapely.Point([x['LONGITUDE(deg)'], x['LATITUDE(deg)']]), axis=1)
    valid_pts.columns = ['longitude','latitude','geometry']



def read_tempo_coords(file_name="tempo_coord_array.pkl", folder_name="./"):

    ## preload TEMPO coordinate matrix (refer to "/download/tempo.py")
    with open(folder_name + file_name, "rb") as f:
        coord = pickle.load(f)
        
    coord_gdf = gpd.GeoDataFrame({
        "longitude": coord[:,:,0].flatten(),
        "latitude": coord[:,:,1].flatten(),
    })



def process_pixel_extraction(file_list,                         ##  TEMPO NO2 file list
                             support_file_list,                 ##  TEMPO support data file list
                             folder,                            ##  TEMPO file folder
                             main_col="Tropo_NO2",              ##  TEMPO Tropospheric NO2 VCD field name
                             qa_flag_col=None,                  ##  specify QA flag field name if applying QA screening 
                             cloud_col='eff_cloud_fraction',    ##  cloud fraction variable if applying QA screening
                             cloud_fraction=0.20,               ##  cloud fraction cutoff if applying QA screening
                             snow_col='snow_ice_fraction',      ##  snow/ice fraction variable if applying QA screening
                             snow_fraction=0.05,                ##  snow/ice fraction cutoff if applying QA screening
                             cutoff_minute=59                   ##  use for temporal colocation, adjust hours accounting for difference between scanning start and actual observations over targeted areas
                            ):
    
    extracted = []
    
    for obs, aux in zip(file_list, support_file_list):
        
        ### 1) Load data
        
        with open(folder + obs, "rb") as f:
            case = pickle.load(f)
            
        with open(folder + aux, "rb") as f:
            case_support = pickle.load(f)
            
        ### 2) Process data
            
        for n in range(len(case)):
            
            tem = copy.deepcopy(case[n][main_col][0,:,:])
            
            if qa_flag_col != None:
                screen = copy.deepcopy(case[n][qa_flag_col][0,:,:])
                # screening with QA flag
                tem[screen != 0] = np.nan
            
            # remove high cloud fraction > 0.20
            cloud_screen = copy.deepcopy(case_support[n][cloud_col][0,:,:])
            tem[cloud_screen > cloud_fraction] = np.nan
            
            # remove high snow/ice fraction > 0.05
            snow_screen = copy.deepcopy(case_support[n][snow_col][0,:,:])
            tem[snow_screen > snow_fraction] = np.nan
            
            # inverting along horizontal axis
            tem = np.flip(tem, axis=0)

        ### 3) Geo-Locate data 
                
            case_df = gpd.GeoDataFrame({
                "longitude": coord_gdf.longitude,
                "latitude": coord_gdf.latitude,
                "observation": tem.flatten()
            })
            
            case_xr = xr.Dataset.from_dataframe(case_df.set_index(['longitude','latitude']))
            
            valid_pts_tempo = []
            
            for p in range(valid_pts.shape[0]):
                
                #  linearly interpolating TEMPO observations to the exact monitoring locations
                valid_pt_tempo_val = case_xr.interp(longitude = valid_pts.longitude.tolist()[p], 
                                                    latitude = valid_pts.latitude.tolist()[p], 
                                                    method = "linear")
                #  add to collection
                valid_pts_tempo.append(valid_pt_tempo_val.observation.values)
            
            valid_pts_at_n = copy.deepcopy(valid_pts)
            valid_pts_at_n['fem_tempo_extract'] = pd.Series(valid_pts_tempo)
            valid_pts_at_n["time"] = str(case[n]["Time"])[2:-2]
            
            extracted.append(valid_pts_at_n)
            
    extracted = pd.concat(extracted)
            
    return extracted



def process_time_groups(extracted_obs):

    ## after running the above function and compiling an extracted dataframe, retrieve temporal groupings
    
    extracted_obs_tempo['time'] = pd.to_datetime(extracted_obs_tempo['time'])
    extracted_obs_tempo['minute'] = extracted_obs_tempo['time'].apply(lambda x: x.minute)
    extracted_obs_tempo['week'] = extracted_obs_tempo['time'].apply(lambda x: x.week)
    extracted_obs_tempo['month'] = extracted_obs_tempo['time'].apply(lambda x: x.month)
    extracted_obs_tempo['weekday'] = extracted_obs_tempo['time'].apply(lambda x: x.weekday())
    extracted_obs_tempo['hour'] = extracted_obs_tempo['time'].apply(lambda x: x.hour)
    extracted_obs_tempo['date'] = extracted_obs_tempo['time'].apply(lambda x: datetime.datetime.strftime(x, "%Y-%m-%d"))
    extracted_obs_tempo['Timestamp(UTC)'] = extracted_obs_tempo.apply(lambda x: x['date'] + "T" + str(x['hour_correct']) + ":00:00-0000", axis=1)



def merge_tempo_airnow(extracted_obs_tempo, valid_monitor_obs):

    #  ensure no duplicated hours for each monitor, otherwise taking the average
    extracted_obs_tempo = extracted_obs_tempo.groupby(
        ['Timestamp(UTC)','hour_correct','date','weekday','week','month','longitude', 'latitude']
    ).agg({'fem_tempo_extract': np.mean}).reset_index()

    #  rounding coordinates for effective merging (avoiding long float leading to mismatch)
    extracted_obs_tempo['longitude'] = extracted_obs_tempo['longitude'].apply(lambda x: round(x, 2))
    extracted_obs_tempo['latitude'] = extracted_obs_tempo['latitude'].apply(lambda x: round(x, 2))
    valid_monitor_obs['LONGITUDE(deg)'] = valid_monitor_obs['LONGITUDE(deg)'].apply(lambda x: round(x, 2))
    valid_monitor_obs['LATITUDE(deg)'] = valid_monitor_obs['LATITUDE(deg)'].apply(lambda x: round(x, 2))

    #  merging
    extracted_obs_tempo = extracted_obs_tempo.merge(valid_monitor_obs[['Timestamp(UTC)', 'LONGITUDE(deg)', 'LATITUDE(deg)','value']],
                                                    left_on = ['Timestamp(UTC)','longitude','latitude'],
                                                    right_on = ['Timestamp(UTC)', 'LONGITUDE(deg)', 'LATITUDE(deg)'],
                                                    how = "left")



def retain_CA_stations(extracted_obs_tempo):

    # read Canadian provincial boundaries
    path = "./Shapefiles/"
    prov_shp = gpd.read_file(path + "Province/lpr_000a21a_e.shp")
    prov_shp = prov_shp[prov_shp["PRUID"].isin(["10","11","12","13","35","24"])]        # 6 provinces in eastern Canada
    prov_shp = prov_shp.to_crs(epsg=4326)

    # remove monitoring stations in the US 
    extracted_obs_tempo['geometry'] = extracted_obs_tempo.apply(lambda x: shapely.Point(x['longitude'], x['latitude']), axis=1)
    extracted_obs_tempo = gpd.GeoDataFrame(extracted_obs_tempo)
    extracted_obs_tempo = extracted_obs_tempo.set_crs(epsg=4326)

    extracted_obs_tempo = gpd.sjoin(left_df = extracted_obs_tempo, right_df = prov_shp[['PRUID','geometry']], how='left') 
    extracted_obs_tempo = extracted_obs_tempo[~pd.isnull(extracted_obs_tempo['PRUID'])].reset_index(drop=True)



def temporal_aggregate(extracted_obs_tempo, time_col):

    extracted_obs_tempo_aggregate = extracted_obs_tempo.groupby([time_col,'longitude','latitude']).agg({
        k: [np.nanmean, np.nanstd] for k in 
            ['fem_tempo_extract','fem_value']
    }).reset_index()

    extracted_obs_tempo_aggregate.columns = [time_col,'longitude','latitude'] + list(
                                        chain(*[
                                                [x, x + "_err"] for x in [
                                                    'fem_tempo_extract','fem_value'
                                                ]
                                            ])
                                    )
    
    return extracted_obs_tempo_aggregate



def cosine_similarity(ts1, ts2):

    """Calculates the cosine similarity between two time series."""
    ts1 = np.array(ts1)
    ts2 = np.array(ts2)
    if len(ts1) != len(ts2):
      raise ValueError("Time series must have the same length for cosine similarity.")
    
    return dot(ts1, ts2) / (norm(ts1) * norm(ts2))



def rmse(ts1, ts2):
    
    """Calculates the Euclidean distance between two time series."""
    ts1 = np.array(ts1)
    ts2 = np.array(ts2)
    if len(ts1) != len(ts2):
        raise ValueError("Time series must have the same length for Euclidean distance.")

    return np.sqrt(np.mean((ts1 - ts2) ** 2))



def plot_regression(extracted_obs_tempo_aggregate):

    plt.rcParams["figure.figsize"] = [12.8, 9.6]

    # add error bar estimates
    plt.errorbar(data = extracted_obs_tempo_aggregate[["fem_tempo_extract","fem_tempo_extract_err","fem_value","fem_value_err"]].dropna(), 
                 x="fem_tempo_extract", y="fem_value", 
                 xerr="fem_tempo_extract_err", yerr="fem_value_err", 
                 elinewidth=0.12, capsize=0, capthick=0, markeredgewidth=0.7, fmt='-',
                 linestyle='', color='grey'
                )

    # add regression lines with scatter points
    seaborn.regplot(data = extracted_obs_tempo_aggregate[["fem_tempo_extract","fem_tempo_extract_err","fem_value","fem_value_err"]].dropna(), 
                    x="fem_tempo_extract", y="fem_value",
                    marker='o', line_kws={"color": "red"}, scatter_kws={'s':20, "color": "black"}
                    )

    plt.xlabel("TEMPO " + "$NO_{2}$" + " VCD (1 x " + "$10^{15}$ " + "$  molec/cm^{2})$", fontsize=16)
    plt.ylabel("Surface " + "$NO_{2}$" + " measurement (ppb)", fontsize=16)
    plt.text(np.nanmax(extracted_obs_tempo_aggregate["fem_tempo_extract"])*0.95, np.nanmax(extracted_obs_tempo_aggregate["fem_value"])*0.95, 
            "Spearman's R = " + str(round(extracted_obs_tempo_aggregate[['fem_tempo_extract','fem_value']].corr(method="spearman").iloc[0,1], 2)),
            fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)



def plot_trend(extracted_obs_tempo_aggregate, time_col, plot_mode):

    #  min-max rescaling the values due to differences in SI units and ranges (molec/cm2 vs ppb)
    extracted_obs_tempo_aggregate['fem_tempo_extract'] = MinMaxScaler().fit_transform(extracted_obs_tempo_aggregate['fem_tempo_extract'].values.reshape(-1,1))
    extracted_obs_tempo_aggregate['fem_value'] = MinMaxScaler().fit_transform(extracted_obs_tempo_aggregate['fem_value'].values.reshape(-1,1))

    #  plotting dataframe
    plot_df = extracted_obs_tempo_aggregate[['fem_tempo_extract','fem_value',"week"]]
    plot_df.columns = ['TEMPO', 'AirNow', time_col.capitalize()]

    #  time series plot
    #  either expressed in lines or boxes & whiskers
    assert plot_mode in ["line", "box"], "avaialble choices are ['line', 'box']"

    plt.rcParams["figure.figsize"] = [6.4, 4.8]

    if plot_mode == "line":

        g = seaborn.lineplot(data = pd.melt(plot_df[['TEMPO', 'AirNow', time_col.capitalize()]], [time_col.capitalize()]), 
                            x = time_col.capitalize(), 
                            y = 'value', 
                            hue = 'variable')

    elif plot_mode == "box":

        g = seaborn.boxplot(data = pd.melt(plot_df [['TEMPO', 'AirNow', time_col.capitalize()]], [time_col.capitalize()]), 
                            x = time_col.capitalize(), 
                            y = 'value', 
                            hue = 'variable',
                            palette="pastel", width=.5)

    plt.ylabel("Normalized Values (Min-Max scaling)")
    plt.legend(title="")

    temp = g.xaxis.get_ticklabels()
    temp = list(set(temp) - set(temp[::4]))
    for label in temp:
        label.set_visible(False)
