import pandas as pd
import geopandas as gpd
import numpy as np
import os
import pickle



'''  Initialization  '''

#  Atlantic Canada
folder = f"./Atlantic Canada/Output/"
no2_agg = sorted([x for x in os.listdir(folder) if "no2_" in x and "agg" in x])
no2_obs_counts = sorted([x for x in os.listdir(folder) if "no2_" in x and "obs_count" in x])

#  QW corridor
folder_qw = f"./QW corridor/Output/"
no2_agg_qw = sorted([x for x in os.listdir(folder_qw) if "no2_" in x and "agg" in x])
no2_obs_counts_qw = sorted([x for x in os.listdir(folder_qw) if "no2_" in x and "obs_count" in x])

# 2023 Sept - 2024 Aug
season_keys = {
        "SON": ["2023_09", "2023_10", "2023_11"],
        "DJF": ["2023_12", "2024_01", "2024_02"],
        "MAM": ["2024_03", "2024_04", "2024_05"],
        "JJA": ["2024_06", "2024_07", "2024_08"]
}

# invalid hours from UTC to time zones in eastern Canada
invalid_hours = ["-4","-3","-2"]
invalid_hours_corrected_atl = ["20","21","22"]
invalid_hours_corrected_qw = ["19","20","21"]




'''  Prerequisite:  refer to "qa_aggregate.py" on how to get temporally aggregated arrays  '''



def load_monthly_mean(no2_agg: list, no2_obs_counts: list):

    ##  1)  Monthly mean NO2

    no2_monthly = {}

    for n in [y for y in no2_agg if "month_agg" in y]:
        # read
        with open(folder + n, "rb") as f:
            agg_data = pickle.load(f)
        # month labels (keys)
        no2_monthly[n[4:11]] = agg_data[next(iter(agg_data))]

    ##  2)  Monthly available observations from QA/QC control

    no2_monthly_counts = {}

    for n in [y for y in no2_obs_counts if "month_obs_count" in y]:
        # read
        with open(folder + n, "rb") as f:
            agg_data_count = pickle.load(f)
        # month labels (keys)
        no2_monthly_counts[n[4:11]] = agg_data_count[next(iter(agg_data_count))]

    return no2_monthly, no2_monthly_counts



def get_annual_mean(no2_monthly: dict, no2_monthly_counts: dict):

    ##  Weighted mean from monthly averages --> annual average

    weighted_mean_numerator = np.nansum([
        no2_monthly[k] * no2_monthly_counts[k] for k in no2_monthly.keys()
    ], 
    axis=0)

    weighted_mean_denominator = np.nansum([
        no2_monthly_counts[k] for k in no2_monthly.keys()
    ], 
    axis=0)

    no2_annual = weighted_mean_numerator / weighted_mean_denominator

    return no2_annaul



def get_seasonal_mean(no2_monthly: dict, no2_monthly_counts: dict, season_keys: dict):

    ##  Weighted mean from monthly averages --> seasaonl average

    no2_seasonal = {}
    no2_seasonal_count = {}

    for k in season_keys.keys():

        tem = [no2_monthly[q] for q in season_keys[k]]
        tem_count = [no2_monthly_counts[q] for q in season_keys[k]]

        season_numerator = np.nansum([obs * count for obs, count in zip(tem, tem_count)], axis=0)
        season_denominator = np.nansum(tem_count, axis=0)

        no2_seasonal[k] = season_numerator / season_denominator
        no2_seasonal_count[k] = season_denominator

    return no2_seasonal, no2_seasonal_count



def load_weekday_mean(no2_agg: list, no2_obs_counts: list):

    ##  1)  Weekday (Monday - Sunday) mean NO2 per each month

    no2_wkday = {}

    for n in [y for y in no2_agg if "weekdays_agg" in y]:
        # read
        with open(folder + n, "rb") as f:
            agg_data = pickle.load(f)
        # weekday labels (keys)
        month_wkday = {}
        for w in agg_data.keys():
            month_wkday[w] = agg_data[w]
        # month labels (keys)
        no2_wkday[n[4:11]] = month_wkday

    ##  2)  Weekly (Monday - Sunday) available observations from QA/QC control per each month

    no2_wkday_counts = {}

    for n in [y for y in no2_obs_counts if "weekdays_obs_count" in y]:
        # read
        with open(folder + n, "rb") as f:
            agg_data_count = pickle.load(f)
        # weekday labels (keys)
        month_wkday_counts = {}
        for w in agg_data_count.keys():
            month_wkday_counts[w] = agg_data_count[w]
        # month labels (keys)
        no2_wkday_counts[n[4:11]] = month_wkday_counts

    return no2_wkday, no2_wkday_counts



def get_monthly_weekday_weekend_mean(no2_wkday: dict, no2_wkday_counts: dict):
    
    no2_wkday_binary = {}
    no2_wkday_binary_counts = {}

    for k in no2_wkday.keys():

        tem = {}
        tem_count = {}

        ##  weighted mean for weekdays (Mon - Fri)

        wkday_weighted_mean_numerator = np.nansum([
            no2_wkday[k][w] * no2_wkday_counts[k][w] for w in [0, 1, 2, 3, 4]
        ], 
        axis=0)

        wkday_weighted_mean_denominator = np.nansum([
            no2_wkday_counts[k][w] for w in [0, 1, 2, 3, 4]
        ], 
        axis=0)

        wkday_agg = weighted_mean_numerator / weighted_mean_denominator

        ##  weighted mean for weekends (Sat - Sun)

        wkend_weighted_mean_numerator = np.nansum([
            no2_wkday[k][w] * no2_wkday_counts[k][w] for w in [5, 6]
        ], 
        axis=0)

        wkend_weighted_mean_denominator = np.nansum([
            no2_wkday_counts[k][w] for w in [5, 6]
        ], 
        axis=0)

        wkend_agg = wkend_weighted_mean_numerator / wkend_weighted_mean_denominator

        ##  write to new dict

        tem["weekday"] = wkday_agg
        tem["weekend"] = wkend_agg
        tem_count["weekday"] = wkday_weighted_mean_denominator
        tem_count["weekend"] = wkend_weighted_mean_denominator

        no2_wkday_binary[k] = tem
        no2_wkday_binary_counts[k] = tem_count
    
    return no2_wkday_binary, no2_wkday_binary_counts



def get_multimonth_weekday_weekend_mean(no2_wkday_binary: dict, no2_wkday_binary_counts: dict):

    ##  Weighted mean from monthly weekday & weekend averages --> annual/seasonal weekday & weekend average

    no2_wkday_binary_multi = {}

    for k in ["weekday", "weekend"]:

        weighted_mean_numerator = np.nansum([
            no2_wkday_binary[q][k] * no2_wkday_binary_counts[q][k] for q in no2_wkday_binary.keys()
        ], 
        axis=0)

        weighted_mean_denominator = np.nansum([
            no2_wkday_binary_counts[q][k] for q in no2_wkday_binary.keys()
        ], 
        axis=0)

        agg = weighted_mean_numerator / weighted_mean_denominator

        no2_wkday_binary_multi[k] = agg

    return no2_wkday_binary_multi



def get_seasonal_weekday_weekend_mean(no2_wkday: dict, no2_wkday_counts: dict, season_keys: dict):

    no2_wkday_seasonal = {}

    for k in season_keys.keys():

        # read data of the three months in the season
        season_tem = {q: no2_wkday[q] for q in season_keys[k]}
        season_tem_counts = {q: no2_wkday_counts[q] for q in season_keys[k]}

        # get monthly mean for the three months in the season
        no2_wkday_binary, no2_wkday_binary_count = get_monthly_weekday_weekend_mean(season_tem, season_tem_counts)

        ##  Weighted mean from monthly weekday & weekend averages --> seasonal weekday & weekend average

        no2_wkday_seasonal[k] = get_multimonth_weekday_weekend_mean(no2_wkday_binary, no2_wkday_binary_count)

    return no2_wkday_seasonal



def get_annual_weekday_weekend_mean(no2_wkday_binary: dict, no2_wkday_binary_count: dict):

    ##  Weighted mean from monthly weekday & weekend averages --> annual weekday & weekend average

    return get_multimonth_weekday_weekend_mean(no2_wkday_binary, no2_wkday_binary_count)



def load_weekday_hourly_mean(no2_agg: list, no2_obs_counts: list):

    no2_dayhour = {}

    for n in [y for y in no2_agg if "weekday_hour_agg" in y]:
        # read
        with open(folder + n, "rb") as f:
            agg_data = pickle.load(f)
        # weekday + hour labels (keys)
        month_dayhour = {}
        for w in agg_data.keys():
            month_dayhour[w] = agg_data[w]
        # month labels (keys)
        no2_dayhour[n[4:11]] = month_dayhour

    no2_dayhour_counts = {}

    for n in [y for y in no2_obs_counts if "weekday_hour_obs_count" in y]:
        # read
        with open(folder + n, "rb") as f:
            agg_data_count = pickle.load(f)
        # weekday + hour labels (keys)
        month_dayhour_counts = {}
        for w in agg_data_count.keys():
            month_dayhour_counts[w] = agg_data_count[w]
        # month labels (keys)
        no2_dayhour_counts[n[4:11]] = month_dayhour_counts

    return no2_dayhour, no2_dayhour_counts



def rearrange_weekday_hours(no2_dayhour: dict, no2_dayhour_counts: dict):

    ##  Splitting the independent "weekday + hour" keys into hierarchical dictionary keys of "weekday" --> "hour" 

    for k in no2_dayhour.keys():
        
        ##  retrieve all hour interval labels from the dictionary keys
        hours = sorted(list(set([str(x.split("_hour")[1]).zfill(2) for x in no2_dayhour[k].keys()])))

        hour_dict = {}
        hour_count_dict = {}

        for j in hours:

            if len([i for i in no2_dayhour[k].keys()
                        if "hour" + str(int(j)) in i and ("day0" in i or "day1" in i or "day2" in i or "day3" in i or "day4" in i)]) > 0 \
                and \
                len([i for i in no2_dayhour[k].keys()
                        if "hour" + str(int(j)) in i and ("day5" in i or "day6" in i)]) > 0:

                bin_dict = {}
                bin_count_dict = {}

                ##  weighted mean for weekdays (Mon - Fri)

                wkday_numerator = np.nansum([
                    no2_dayhour[k][i] * no2_dayhour_counts[k][i]
                    for i in no2_dayhour[k].keys()
                    if "hour" + str(int(j)) in i and ("day0" in i or "day1" in i or "day2" in i or "day3" in i or "day4" in i)
                ], axis=0)

                wkday_denominator = np.nansum([
                    no2_dayhour_counts[k][i]
                    for i in no2_dayhour[k].keys()
                    if "hour" + str(int(j)) in i and ("day0" in i or "day1" in i or "day2" in i or "day3" in i or "day4" in i)
                ], axis=0)

                bin_dict["weekday"] = wkday_numerator / wkday_denominator
                bin_count_dict["weekday"] = wkday_denominator

                ##  weighted mean for weekends (Sat - Sun)

                wkend_numerator = np.nansum([
                    no2_dayhour[k][i] * no2_dayhour_counts[k][i]
                    for i in no2_dayhour[k].keys()
                    if "hour" + str(int(j)) in i and ("day5" in i or "day6" in i)
                ], axis=0)

                wkend_denominator = np.nansum([
                    no2_dayhour_counts[k][i]
                    for i in no2_dayhour[k].keys()
                    if "hour" + str(int(j)) in i and ("day5" in i or "day6" in i)
                ], axis=0)

                bin_dict["weekend"] = wkend_numerator / wkend_denominator
                bin_count_dict["weekend"] = wkend_denominator

                ##  write to new dictionary

                hour_dict[j] = bin_dict
                hour_count_dict[j] = bin_count_dict

        no2_dayhour[k] = hour_dict
        no2_dayhour_counts[k] = hour_count_dict



def clean_invalid_hours(no2_dayhour: dict, no2_dayhour_counts: dict,
                        invalid_hours: list, invalid_hours_corrected: list):

    for j in no2_dayhour.keys():

        ##  Due to raw timestamps in UTC, negative hour labels may be obtained when adjusting time zones
        ##  Combine/Add those to the observations of the correct hour labels

        for k in zip(invalid_hours, invalid_hours_corrected):

            if k[0] in no2_dayhour[j].keys():

                if k[1] in no2_dayhour[j].keys():

                    for w in ["weekday", "weekend"]:

                        numerator = np.nansum([
                            no2_dayhour[j][k[0]][w] * no2_dayhour_counts[j][k[0]][w],
                            no2_dayhour[j][k[1]][w] * no2_dayhour_counts[j][k[1]][w]
                        ], 
                        axis=0)

                        denominator = np.nansum([
                            no2_dayhour[j][k[0]][w],
                            no2_dayhour[j][k[1]][w]
                        ], 
                        axis=0)

                        no2_dayhour[j][k[1]][w] = numerator / denominator
                        no2_dayhour_counts[j][k[1]][w] = denominator

                    del no2_dayhour[j][k[0]]
                    del no2_dayhour_counts[j][k[0]]

                else:

                    no2_dayhour[j][k[1]] = no2_dayhour[j][k[0]]
                    no2_dayhour_counts[j][k[1]] = no2_dayhour_counts[j][k[0]]

                    del no2_dayhour_progress[j][k[0]]
                    del no2_dayhour_counts[j][k[0]]



def get_weekday_weekend_hourly_mean(no2_dayhour: dict, no2_dayhour_counts: dict):

    no2_dayhour_agg = {}
    no2_dayhour_count_agg = {}

    hours = set.union(*map(set, [list(no2_dayhour[k].keys()) for k in no2_dayhour.keys()]))

    for h in hours:

        no2_dayhour_agg[h] = {}
        no2_dayhour_count_agg[h] = {}

        ##  weighted mean for weekdays (Mon - Fri)

        wkday_numerator = np.nansum([
            no2_dayhour[k][h]["weekday"] * no2_dayhour_counts[k][h]["weekday"] 
            for k in no2_dayhour.keys() 
            if h in no2_dayhour[k].keys()
        ], 
        axis=0)

        wkday_denominator = np.nansum([
            no2_dayhour_counts[k][h]["weekday"] 
            for k in no2_dayhour.keys() 
            if h in no2_dayhour[k].keys()
        ], 
        axis=0)

        agg_wkday = wkday_numerator / wkday_denominator
        count_wkday = wkday_denominator

        ##  weighted mean for weekends (Sat - Sun)

        wkend_numerator = np.nansum([
            no2_dayhour[k][h]["weekend"] * no2_dayhour_counts[k][h]["weekend"] 
            for k in no2_dayhour.keys() 
            if h in no2_dayhour[k].keys()
        ], 
        axis=0)

        wkend_denominator = np.nansum([
            no2_dayhour_counts[k][h]["weekend"] 
            for k in no2_dayhour.keys() 
            if h in no2_dayhour[k].keys()
        ], 
        axis=0)

        agg_wkend = wkend_numerator / wkend_denominator
        count_wkend = wkend_denominator

        ##  write to new dictionary

        no2_dayhour_agg[h]['weekday'] = agg_wkday
        no2_dayhour_agg[h]['weekend'] = agg_wkend

        no2_dayhour_count_agg[h]['weekday'] = count_wkday
        no2_dayhour_count_agg[h]['weekend'] = count_wkend

    return no2_dayhour_agg, no2_dayhour_count_agg



def get_seasonal_weekday_weekend_hourly_mean(no2_dayhour: dict, no2_dayhour_counts: dict, season_keys: dict,
                                             invalid_hours: list, invalid_hours_corrected: list):

    ##  Weighted mean from monthly weekday & weekend hourly averages --> seasonal weekday & weekend hourly average

    no2_seasonal_dayhour = {}
    no2_seasonal_dayhour_count = {}

    for s in season_keys.keys():

        no2_dayhour = {q: no2_dayhour[q] for q in season_keys[s]}
        no2_dayhour_count = {q: no2_dayhour_counts[q] for q in season_keys[s]}

        rearrange_weekday_hours(no2_dayhour, no2_dayhour_counts)
        clean_invalid_hours(no2_dayhour, no2_dayhour_counts, invalid_hours, invalid_hours_corrected)

        no2_dayhour_agg, no2_dayhour_count_agg = get_weekday_weekend_hourly_mean(no2_dayhour, no2_dayhour_counts)

        no2_seasonal_dayhour[s] = no2_dayhour_agg
        no2_seasonal_dayhour_count[s] = no2_dayhour_count_agg

    return no2_seasonal_dayhour, no2_seasonal_dayhour_count
