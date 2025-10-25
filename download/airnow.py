import pandas as pd
import time
import pyrsig
from datetime import datetime


#  attribute name
features_set = [
    'airnow.no2'
]


#  batches / divisions of the study period
download_dates = [
    pd.date_range(start="2023-09-01", end="2024-08-31").to_pydatetime(),
]


##  Atlantic Canada
atlantic_bbox = (-69.05, 43.05, -52.05, 51.05)

##  Quebec City - Windsor corridor
qw_bbox = (-86.05, 41.20, -69.05, 48.20)



## Downloading function for execution

def downloading(region_bbox):

    airNow_payload = []

    for dates in download_dates:

        for date in dates:

            rsigapi = pyrsig.RsigApi(
                bdate = datetime.strftime(date, "%Y-%m-%d"),
                bbox = region_bbox
                )

            sub_airNow_payload = []

            for f in features_set:
                df = rsigapi.to_dataframe(f)
                df['pollutant'] = df.columns[-2]
                df.columns = list(df.columns[0:-3]) + ["value"] + list(df.columns[-2:])

                sub_airNow_payload.append(df)

            airNow_payload.append(pd.concat(sub_airNow_payload))

        pd.concat(sub_airNow_payload).to_excel(f"/AirNow/NO2/" + datetime.strftime(date, "%Y-%m-%d") + ".xlsx", index=False)

        if len(airNow_payload) > 0:
            print("completed run for " + datetime.strftime(date, "%Y-%m-%d"))

        time.sleep(60)

