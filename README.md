###  ETL and visualization workflow for spatiotemporal analysis of TEMPO L3 Tropospheric NO<sub>2</sub> in Ontario, Quebec, and Atlantic Canada
---------------------------------------------------------------------------------------------------------

This repository presents the Python codes used for the analysis of the manuscript, <i>"Tropospheric NO<sub>2</sub> Patterns in Eastern Canada Using the First Year of TEMPO Observations"</i>, submitted to JGR: Atmosphere. 
All TEMPO data were sourced from earthdata.nasa.gov (by NASA Atmospheric Science Data Center), provisional gridded Level 3 tropospheric NO<sub>2</sub> column (DOI: <a her="https://doi.org/10.5067/IS-40e/TEMPO/NO2_L3.003">https://doi.org/10.5067/IS-40e/TEMPO/NO2_L3.003</a>) .

Refer to the following scripts in this repository and their respective purposes:

- ```/download/```
  - ```tempo.py``` : subsetting and extraction of required attributes from TEMPO netCDF files
  - ```tropomi.py``` : downloading process of TROPOMI Level 3 data from Google Earth Engine
  - ```airnow.py``` : downloading process of hourly measurements of surface monitors from AirNow via RSIG

- ```qa_aggregate.py``` : masking of pixel observations with high uncertainties and temporally aggregated data arrays
- ```preload_aggregate.py``` : fast-loading the temporally aggregated data arrays for visualization

- ```diurnal.ipynb``` : plots assessing city-level diurnal trends with stratification by seasons and weekdays/weekends
- ```monitoring_gaps.ipynb``` : plots assessing under-monitored communities based on TEMPO spatial variability
- ```columnar_validate.py``` :  spatiotemporally colocated TEMPO and TROPOMI comparison
- ```column_surface.py``` :  spatiotemporally colocated TEMPO & AirNow observations at surface monitor locations
