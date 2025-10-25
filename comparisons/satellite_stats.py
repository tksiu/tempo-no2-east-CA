from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
import seaborn
import pickle
import copy
import numpy as np
import pandas as pd
import xarray as xr
import scipy
import statsmodels.api as sm



#  colocated and projected at 2 km x 2 km resolutions (TEMPO grids)

from satellite_colocation import spatial_matching

atl_colocation_instance = spatial_matching(tempo_coord_folder = "./coordinates/", 
                                            tempo_coord_file = "coord_array_atl.pkl",
                                            tropomi_coord_folder = "./coordinates/",
                                            tropomi_coord_file = "coord_array_tropomi_atl.pkl")
qw_colocation_instance = spatial_matching(tempo_coord_folder = "./coordinates/", 
                                            tempo_coord_file = "coord_array_qw.pkl",
                                            tropomi_coord_folder = "./coordinates/",
                                            tropomi_coord_file = "coord_array_tropomi_qw.pkl")

coord_atl_gdf = atl_colocation_instance.tempo_coords_gdf
coord_atl = atl_colocation_instance.tempo_coords_xr

coord_qw_gdf = qw_colocation_instance.tempo_coords_gdf
coord_qw = qw_colocation_instance.tempo_coords_xr



#  read intermediate processed datasets

with open("./columnar_spatiotemporally_colocated/atl_tempo.pkl", "rb") as f:
    tempo_data_atl = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/atl_tropomi.pkl", "rb") as f:
    tropomi_data_atl = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/atl_tempo_timestamp.pkl", "rb") as f:
    tempo_match_time_atl = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/atl_tropomi_timestamp.pkl", "rb") as f:
    tropomi_match_time_atl = pickle.load(f)

with open("./columnar_spatiotemporally_colocated/atl_terrain_height.pkl", "rb") as f:
    terrain_height_atl = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/atl_albedo.pkl", "rb") as f:
    surface_albedo_atl = pickle.load(f)

with open("./columnar_spatiotemporally_colocated/qw_tempo.pkl", "rb") as f:
    tempo_data_qw = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/qw_tropomi.pkl", "rb") as f:
    tropomi_data_qw = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/qw_tempo_timestamp.pkl", "rb") as f:
    tempo_match_time_qw = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/qw_tropomi_timestamp.pkl", "rb") as f:
    tropomi_match_time_qw = pickle.load(f)

with open("./columnar_spatiotemporally_colocated/qw_terrain_height.pkl", "rb") as f:
    terrain_height_qw = pickle.load(f)
    
with open("./columnar_spatiotemporally_colocated/qw_albedo.pkl", "rb") as f:
    surface_albedo_qw = pickle.load(f)



#  accommodate numpy arrays to xarray datasets

def fetch_plot_Xr(aggregate, coords):

    plot_gdf = copy.deepcopy(coords)
    plot_gdf["observation"] = aggregate.flatten()

    plot_xr = xr.Dataset.from_dataframe(plot_gdf.set_index(['longitude','latitude']))
    plot_xr = plot_xr.rio.write_crs("epsg:4326")
    plot_xr = plot_xr.rename({"longitude": "x", "latitude": "y"})
    plot_xr = plot_xr.transpose("y", "x")
    plot_xr = plot_xr.reindex(y=list(reversed(plot_xr.y)))

    return plot_xr

# remove observations < 1e14 (extremely low values, only less than 1% of data)
for x in range(len(tempo_data_atl)):
    tempo_data_atl[x][tempo_data_atl[x] < 1e14] = np.nan
    tropomi_data_atl[x][tropomi_data_atl[x] < 1e14] = np.nan

ny_exclusion = 45     ## exlcude New York city from the spatial extent
for t in range(len(tempo_data_qw)):
    tempo_data_qw[t] = tempo_data_qw[t][:-ny_exclusion,:]

for x in range(len(tempo_data_qw)):
    tempo_data_qw[x][tempo_data_qw[x] < 1e14] = np.nan
    tropomi_data_qw[x][tropomi_data_qw[x] < 1e14] = np.nan

tempo_agg = fetch_plot_Xr(np.nanmean(tempo_data_atl, axis=0), coord_atl_gdf)
tropomi_agg = fetch_plot_Xr(np.nanmean(tropomi_data_atl, axis=0), coord_atl_gdf)

tempo_agg_qw = fetch_plot_Xr(np.nanmean(tempo_data_qw, axis=0), coord_qw_gdf)
tropomi_agg_qw = fetch_plot_Xr(np.nanmean(tropomi_data_qw, axis=0), coord_qw_gdf)



###  Analysis  ###

# add MODIS land cover classes (refer to "/download/modis_lc.py")
lc_atl = pd.read_csv("./land_cover/TEMPO NO2 spatial extents land covers - Atl_reclassified.csv")
lc_qw = pd.read_csv("./land_cover/TEMPO NO2 spatial extents land covers - QW_reclassified.csv")

# gather TEMPO and TROPOMI values for correlation analysis, stratified by land covers
def columnar_analysis_dataframe(tempo_agg, tropomi_agg, land_cover, terrain_height, surface_albedo):
    
    df = pd.DataFrame({
        "tempo": tempo_agg.flatten(),
        "tropomi": tropomi_agg.flatten(),
        "lc": land_cover['reclass'].values,
        "urban_perc": land_cover['urban'].values,
        "water_perc": land_cover['water'].values,
        "vegetation_perc": land_cover['vegetation'].values,
        "terrain_height": np.nanmean(terrain_height, axis=0).flatten(),
        "surface_albedo": np.nanmean(surface_albedo, axis=0).flatten(),
    })
    df = df.dropna()

    df['dev'] = df['tempo'] - df['tropomi']

    df['urban_0.3'] = df['urban_perc'].apply(lambda x: 1 if x > 0.3 else 0)
    df['water_0.3'] = df['water_perc'].apply(lambda x: 1 if x > 0.3 else 0)
    df['vegetation_0.3'] = df['vegetation_perc'].apply(lambda x: 1 if x > 0.3 else 0)

    return df

# gather TEMPO and TROPOMI pixelwise difference mapping
def pixelwise_difference_xr(tempo_agg, tropomi_agg, tempo_coords_xr, percentage=False):

    tempo_coords_xr['tempo_agg'] = (('y', 'x'), tempo_agg)
    tempo_coords_xr['tropomi_agg'] = (('y', 'x'), tropomi_agg)

    tempo_coords_xr = tempo_coords_xr.interpolate_na(dim="y", method="linear").interpolate_na(dim="x", method="linear")

    if percentage:
        diff = (tempo_coords_xr.tempo_agg - tempo_coords_xr.tropomi_agg) / tempo_coords_xr.tropomi_agg * 100
    else:
        diff = tempo_coords_xr.tempo_agg - tempo_coords_xr.tropomi_agg

    return diff

###  Visualization  ###

def columnar_lr_line(columnar_analysis_df):

    # simple linear regression analysis
    sm.OLS(columnar_analysis_df['tropomi'], sm.add_constant(columnar_analysis_df['tempo']))\
        .fit()\
        .summary()

    # 1) scatter plot with density and linear regression line
    seaborn.regplot(data=columnar_analysis_df, x="tempo", y="tropomi",
                    marker='o', line_kws={"color": "red"}, scatter_kws={'s':2, "color": "black"})
    seaborn.kdeplot(data=columnar_analysis_df, x="tempo", y="tropomi", 
                    fill=True)
    plt.xlabel("TEMPO", fontsize=16)
    plt.ylabel("TROPOMI", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.text(np.nanmax(columnar_analysis_df["tempo"])*1.05, np.nanmax(columnar_analysis_df["tropomi"])*1.05, 
            "Pearson's R = {:.2f}".format(scipy.stats.pearsonr(columnar_analysis_df.tempo, columnar_analysis_df.tropomi)[0]),
            "Spearman's R = {:.2f}".format(scipy.stats.spearmanr(columnar_analysis_df.tempo, columnar_analysis_df.tropomi)[0]),
            fontsize=16)

def columnar_lr_residuals(columnar_analysis_df):

    # 2) residual plot of linear regression
    seaborn.residplot(data=columnar_analysis_df, x="tempo", y="tropomi", 
                      scatter_kws={'s':2, "color": "black"})
    plt.xlabel("TEMPO", fontsize=16)
    plt.ylabel("Residual", fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)

def pixel_difference_map(columnar_difference_xr):

    # 3) pixel map for spatial distribution of TEMPO - TROPOMI deviation
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 9))

    southend = 43.05
    northend = 51.05
    westend = -69.05
    eastend = -52.05

    # create a base map
    m = Basemap(projection='merc',
                llcrnrlon = westend,
                llcrnrlat = southend,
                urcrnrlon = eastend,
                urcrnrlat = northend,
                resolution='h',
                ax = ax)
    m.drawcoastlines(color='#808080', linewidth=0.75)
    m.drawcountries(color='#808080', linewidth=1, linestyle=(0, (5, 10)))

    p = m.imshow(columnar_difference_xr.values, origin="upper", cmap='Spectral_r', alpha=0.9, vmax=1.0e15, vmin=-1.0e15)

    plt.title("Difference (TEMPO - TROPOMI) aggregated from the matched images (2 km x 2 km)")
    fig.colorbar(p, ax=ax, location='right', anchor=(0, 0.3), shrink=0.7)


def columnar_bias_regression_on_surface_conditions(columnar_analysis_df):

    # 4) linear regression model investigating how TEMPO - TROPOMI deviation varied by land surface characteristics
    model = sm.OLS(columnar_analysis_df[['dev']] / 1e15, 
                   sm.add_constant(columnar_analysis_df[["urban_perc","water_perc","vegetation_perc","terrain_height","surface_albedo"]])
            ).fit()
    print(model.summary())

    model = sm.OLS(columnar_analysis_df[['dev']] / 1e15, 
                   sm.add_constant(columnar_analysis_df[["urban_0.3","water_0.3","vegetation_0.3","terrain_height","surface_albedo"]])
            ).fit()
    print(model.summary())
