import ee
ee.Authenticate()
ee.Initialize(project='')

import time



## "filterBounds" did not work for Sentinel-5P TROPOMI because each TROPOMI product will cover a global extent
##  define a new function for checking proportions of non-null values (i.e., scene coverage)

def image_perc_mask(image, scale, aoi):
    totPixels = ee.Number(image.unmask(1).reduceRegion(reducer = ee.Reducer.count(),
                                                       scale = scale,
                                                       geometry = aoi,
                                                       ).values().get(0))
    actPixels = ee.Number(image.reduceRegion(reducer = ee.Reducer.count(),
                                             scale = scale,
                                             geometry = aoi,
                                             ).values().get(0))
    percCover = actPixels.divide(totPixels).multiply(100).round()
    return image.set('percCover', percCover)



##  Atlantic Canada

uplon = -52.05
uplat = 51.05
lowlon = -69.05
lowlat = 43.05

bound = [lowlon, lowlat, uplon, uplat]
atl_roi = ee.Geometry.Rectangle(bound)


##  Quebec City - Windsor corridor

uplon = -69.05
uplat = 48.20
lowlon = -86.05
lowlat = 41.20

bound = [lowlon, lowlat, uplon, uplat]
qw_roi = ee.Geometry.Rectangle(bound)



##  single observations (2023 Sept - 2024 Aug)

def download_tropomi_images(roi):

    tropomi_no2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2') \
                    .filterDate('2023-09-01', '2024-09-01') \
                    .filterBounds(roi) \
                    .select(['tropospheric_NO2_column_number_density'])

    tropomi_no2 = tropomi_no2.map(
        lambda image: image_perc_mask(image, 1113.2, roi)).filter(ee.Filter.gte('percCover', 50)
    )

    tropomi_no2_list = tropomi_no2.toList(tropomi_no2.size())

    counter = 0

    while counter >= 0:

        try:
            img = ee.Image(tropomi_no2_list.get(counter))
            date = ee.Date(img.get('system:time_start')).format(None, 'UTC').getInfo()
        except:
            print("End.")
        break

        task = ee.batch.Export.image.toDrive(**{
                'image': img,
                'description': 'TROPOMI_NO2_' + date,
                'folder': "TROPOMI",
                'region': roi,
                'scale': 1113.2,
                'crs':'EPSG:4326',
                'fileFormat':'GeoTIFF',
            })
        task.start() 
        
        if counter == 0:

            ##  coordinate matrix
            coordinates = img.pixelLonLat().select(['longitude','latitude'])

            task = ee.batch.Export.image.toDrive(**{
                    'image': coordinates,
                    'description': 'TROPOMI_NO2_coordinates',
                    'folder': "TROPOMI",
                    'region': roi,
                    'scale': 1113.2,
                    'crs':'EPSG:4326',
                    'fileFormat':'GeoTIFF',
                })
            task.start() 

        time.sleep(1)

        counter += 1
        if counter % 10 == 0:
            print("completed = " + str(counter))



##  1-year Mean (2023 Sept - 2024 Aug)

def download_tropomi_1yr_mean(roi):

    collection = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")\
                            .select("tropospheric_NO2_column_number_density")\
                            .filterDate("2023-09-01", "2024-09-01")

    collection_non_null = collection.map(
        lambda image: image_perc_mask(image, 1113.2, roi)).filter(ee.Filter.gt('percCover', 0)
    )

    mean_tif = collection_non_null.reduce(
        ee.Reducer.mean()
    )

    task = ee.batch.Export.image.toDrive(
        image = mean_tif,
        description = 'Tropomi_NO2_1yr_mean',
        folder ="TROPOMI 2023/09-2024/08",
        region = roi,
        scale = 1113.2,
        crs = 'EPSG:4326',
        fileFormat = 'GeoTIFF',
    )
    task.start()



##  5-year Mean (2018 Sept - 2023 Aug)

def download_tropomi_5yr_mean(roi):

    collection = ee.ImageCollection("COPERNICUS/S5P/OFFL/L3_NO2")\
                            .select("tropospheric_NO2_column_number_density")\
                            .filterDate("2023-09-01", "2024-09-01")

    collection_non_null = collection.map(
        lambda image: image_perc_mask(image, 1113.2, roi)).filter(ee.Filter.gt('percCover', 0)
    )

    mean_tif = collection_non_null.reduce(
        ee.Reducer.mean()
    )

    task = ee.batch.Export.image.toDrive(
        image = mean_tif,
        description = 'Tropomi_NO2_5yr_mean',
        folder ="TROPOMI 2018/09-2023/08",
        region = roi,
        scale = 1113.2,
        crs = 'EPSG:4326',
        fileFormat = 'GeoTIFF',
    )
    task.start()

