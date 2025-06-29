import numpy as np
import pandas as pd
import os
import copy
import pickle
import datetime



###  Extract temporal grouping labels from the timestamp (scanning start time) of the observations

def timestamp_marks(data, 
                    utc_diff_daylight_saving,       ## relative to UTC; e.g., 3 for Atlantic UTC-3 (daylight saving time from March to November)
                    utc_diff_standard,              ## relative to UTC; e.g., 4 for Atlantic UTC-4 (stadndard time from November to March)
                    clock_change_fall,              ## clock change date in fall (November)
                    clock_change_spring,            ## clock change date in spring (March)
                    cutoff_minute=None              ## adjust for lead time between scanning start and actual observations over targeted areas (~ 6-7 min for each granule)
                    ):
    
    timestamp = [pd.to_datetime(x["Time"][0]) for x in data]
    
    # month
    timestamp_month = [x.month for x in timestamp]
    # season
    timestamp_season = ["DJF" if x in [12, 1, 2] else 
                        "MAM" if x in [3, 4, 5] else 
                        "JJA" if x in [6, 7, 8] else 
                        "SON" 
                        for x in timestamp_month]
    # week of year
    timestamp_yearweek = [x.week for x in timestamp]
    # weekday
    timestamp_day = [x.weekday() for x in timestamp]
    # hour (adjusted to Atlantic Time UTC-3/-4 and QW corridor Time UTC-4/-5)
    if cutoff_minute != None:
        timestamp_hour_correct = [x.hour if x.minute < cutoff_minute else 
                                  x.hour + 1 
                                  for x in timestamp]
    else:
        timestamp_hour_correct = [x.hour 
                                  for x in timestamp]
    timestamp_hour_adjust = [timestamp_hour_correct[x] - utc_diff_daylight_saving if timestamp[x] < clock_change_fall or timestamp[x] >= clock_change_spring else 
                             timestamp_hour_correct[x] - utc_diff_standard
                             for x in range(len(timestamp))]
    # secondary dimension
    timestamp_weekday_hour = ["day" + str(timestamp_day[y]) + "_" + "hour" + str(timestamp_hour_adjust[y]) 
                              for y in range(len(timestamp_day))]
    timestamp_week_hour = ["week" + str(timestamp_yearweek[y]) + "_" + "hour" + str(timestamp_hour_adjust[y]) 
                              for y in range(len(timestamp_yearweek))]
    
    return [timestamp,
            timestamp_month, 
            timestamp_season,
            timestamp_yearweek,
            timestamp_day,
            timestamp_hour_adjust,
            timestamp_weekday_hour,
            timestamp_week_hour]



###  Function for extracting one month's data

def temporal_aggregates(timestamp, data, main_col):

    agg = {}
    count = {}
    
    for n in [x for i, x in enumerate(timestamp) if i == timestamp.index(x)]:
        
        ###  fetch all images falling in the specified time groups, e.g., months, weeks, hours
        idx = [x for x in range(len(timestamp)) if timestamp[x] == n]
        fetch = [data[x][main_col].reshape(1, data[x][main_col].shape[0], data[x][main_col].shape[1]) for x in idx]
        fetch_stack = np.concatenate(fetch)
        
        ###  pixelwise mean with non-null values
        fetch_agg = np.nanmean(fetch_stack, axis=0)
        agg[n] = fetch_agg
        
        ###  pixelwise number of images that are valid/good retrieval (not null)
        fetch_valid = np.count_nonzero(~np.isnan(fetch_stack), axis=0)
        count[n] = fetch_valid
    
    return agg, count



###  Function for extracting all data

def process_temporal_aggregates(input_folder,                   # folder storing and reading inputs
                                output_folder,                  # folder for writing aggregated outputs
                                file_list,                      # list of all files in the input folder
                                main_col,                       # target field, i.e., tropospheric NO2 column density
                                use_support_data=False,         # true/false if using support data for QA/QC screening
                                support_file_list=None,         #   if true, specifying list of all support data files
                                qa_flag_col=None,               # main QA flag column
                                cloud_filter=False,             # true/false if applying cloud fraction filtering 
                                cloud_col=None,                 #   if true, specifying the cloud fraction field name
                                cloud_fraction=0.20,            #   if true, specifying the cloud fraction threshold
                                snow_filter=False,              # true/false if applying snow/ice fraction filtering 
                                snow_col=None,                  #   if true, specifying the snow/ice fraction field name
                                snow_fraction=0.05,             #   if true, specifying the snow/ice fraction threshold
                                tempo_minute_cutoff=None,       # whether applying cutoff minutes (lead time between scanning start and actual observations over targeted areas)
                                ):

    for obs, aux in zip(file_list, support_file_list):
        
        ### 1) Load data
        
        with open(input_folder + obs, "rb") as f:
            case = pickle.load(f)
            
        with open(input_folder + aux, "rb") as f:
            case_support = pickle.load(f)
            
        ### 2) Process timestamps
        
        # clock change dates in 2023 - 2024
        adt2ast = datetime.datetime(2023, 11, 5, 0, 0)
        ast2adt = datetime.datetime(2024, 3, 10, 0, 0)

        # get temporal grouping labels
        if tempo_minute_cutoff != None:
            labels =  timestamp_marks(case, adt2ast, ast2adt, tempo_minute_cutoff)
        else:
            labels =  timestamp_marks(case, adt2ast, ast2adt)
        
        timestamp = labels[0]
        timestamp_month = labels[1] 
        timestamp_season = labels[2]
        timestamp_yearweek = labels[3]
        timestamp_day = labels[4]
        timestamp_hour_adjust = labels[5]
        timestamp_weekday_hour = labels[6]
        timestamp_week_hour  = labels[7]
        
        ##  3) Process each array with QA flag screening
        
        for n in range(len(case)):
            
            # read in the data
            tem = copy.deepcopy(case[n][main_col][0,:,:])
            
            if qa_flag_col != None:
                screen = copy.deepcopy(case[n][qa_flag_col][0,:,:])
                # screening with QA flag
                tem[screen != 0] = np.nan
            
            # remove high cloud fraction > 0.2
            if cloud_filter:
                cloud_screen = copy.deepcopy(case_support[n][cloud_col][0,:,:])
                tem[cloud_screen > cloud_fraction] = np.nan
            
            # remove high snow/ice fraction > 0.05
            if snow_filter:
                snow_screen = copy.deepcopy(case_support[n][snow_col][0,:,:])
                tem[snow_screen > snow_fraction] = np.nan
            
            # inverting along horizontal axis
            tem = np.flip(tem, axis=0)
            
            case[n][main_col] = tem
        
        
        month_agg, month_count = temporal_aggregates(timestamp_month, case, main_col)
        yearweek_agg, yearweek_count = temporal_aggregates(timestamp_yearweek, case, main_col)
        day_agg, day_count = temporal_aggregates(timestamp_day, case, main_col)
        hour_agg, hour_count = temporal_aggregates(timestamp_hour_adjust, case, main_col)
        weekday_hour_agg, weekday_hour_count = temporal_aggregates(timestamp_weekday_hour, case, main_col)
        week_hour_agg, week_hour_count = temporal_aggregates(timestamp_week_hour, case, main_col)
        
        
        with open(output_folder + obs.replace(".pkl", "") + "_month_agg.pkl", "wb") as f:
            pickle.dump(month_agg, f)
        with open(output_folder + obs.replace(".pkl", "") + "_month_obs_count.pkl", "wb") as f:
            pickle.dump(month_count, f)
            
        with open(output_folder + obs.replace(".pkl", "") + "_week_agg.pkl", "wb") as f:
            pickle.dump(yearweek_agg, f)
        with open(output_folder + obs.replace(".pkl", "") + "_week_obs_count.pkl", "wb") as f:
            pickle.dump(yearweek_count, f)
    
        with open(output_folder + obs.replace(".pkl", "") + "_weekdays_agg.pkl", "wb") as f:
            pickle.dump(day_agg, f)
        with open(output_folder + obs.replace(".pkl", "") + "_weekdays_obs_count.pkl", "wb") as f:
            pickle.dump(day_count, f)
        
        with open(output_folder + obs.replace(".pkl", "") + "_hour_agg.pkl", "wb") as f:
            pickle.dump(hour_agg, f)
        with open(output_folder + obs.replace(".pkl", "") + "_hour_obs_count.pkl", "wb") as f:
            pickle.dump(hour_count, f)

        with open(output_folder + obs.replace(".pkl", "") + "_weekday_hour_agg.pkl", "wb") as f:
            pickle.dump(weekday_hour_agg, f)
        with open(output_folder + obs.replace(".pkl", "") + "_weekday_hour_obs_count.pkl", "wb") as f:
            pickle.dump(weekday_hour_count, f)
            
        with open(output_folder + obs.replace(".pkl", "") + "_week_hour_agg.pkl", "wb") as f:
            pickle.dump(week_hour_agg, f)
        with open(output_folder + obs.replace(".pkl", "") + "_week_hour_obs_count.pkl", "wb") as f:
            pickle.dump(week_hour_count, f)
        
        print("Processed " + obs.replace(".pkl", ""))



if __name__ == 'main':

    
    ##  Example 1:  Atlantic Canada

    primary_folder = "./Atlantic Canada/"

    input_folders = [primary_folder + "Source/SON/",
                    primary_folder + "Source/DJF/",
                    primary_folder + "Source/MAM/",
                    primary_folder + "Source/JJA/"]

    output_folders = [primary_folder + "Output/SON/",
                    primary_folder + "Output/DJF/",
                    primary_folder + "Output/MAM/",
                    primary_folder + "Output/JJA/"]

    for i, o in zip(input_folders, output_folders):

        no2_obs = [x for x in os.listdir(i) if "_support_data" not in x and "no2" in x and "agg" not in x and "count" not in x]
        no2_support = [x for x in os.listdir(i) if "_support_data" in x and "no2" in x and "agg" not in x and "count" not in x]
        
        process_temporal_aggregates(input_folder = i, 
                                    output_folder = o, 
                                    file_list = no2_obs, 
                                    main_col = 'Tropo_NO2', 
                                    use_support_data = True, 
                                    support_file_list = no2_support,
                                    qa_flag_col = 'main_QA_flag',
                                    cloud_filter=True, 
                                    cloud_col='eff_cloud_fraction', 
                                    cloud_fraction=0.20,
                                    snow_filter=True, 
                                    snow_col='snow_ice_fraction', 
                                    snow_fraction=0.05,
                                    tempo_minute_cutoff = 47, 
                                    negative_values=True)


    ##  Example 2:  QW Corridor

    primary_folder = "./QW corridor/"

    input_folders = [primary_folder + "Source/SON/",
                    primary_folder + "Source/DJF/",
                    primary_folder + "Source/MAM/",
                    primary_folder + "Source/JJA/"]

    output_folders = [primary_folder + "Output/SON/",
                    primary_folder + "Output/DJF/",
                    primary_folder + "Output/MAM/",
                    primary_folder + "Output/JJA/"]

    for i, o in zip(input_folders, output_folders):

        no2_obs = [x for x in os.listdir(i) if "_support_data" not in x and "no2" in x and "agg" not in x and "count" not in x]
        no2_support = [x for x in os.listdir(i) if "_support_data" in x and "no2" in x and "agg" not in x and "count" not in x]
        
        process_temporal_aggregates(input_folder = i, 
                                    output_folder = o, 
                                    file_list = no2_obs, 
                                    main_col = 'Tropo_NO2', 
                                    use_support_data = True, 
                                    support_file_list = no2_support,
                                    qa_flag_col = 'main_QA_flag',
                                    cloud_filter=True, 
                                    cloud_col='eff_cloud_fraction', 
                                    cloud_fraction=0.20,
                                    snow_filter=True, 
                                    snow_col='snow_ice_fraction', 
                                    snow_fraction=0.05,
                                    tempo_minute_cutoff = 34, 
                                    negative_values=True)
