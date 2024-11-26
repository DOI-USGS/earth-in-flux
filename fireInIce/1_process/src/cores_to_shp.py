import pandas as pd
import geopandas as gpd
def main(core_data,core_shpfile):
    cd = pd.read_excel(core_data)
    cd['Latitude'][2] = "58° 51’ 16.1” N" #datafix 2016, C3
    cd['Longitude'][2] = "134° 10’ 31.1” W" #datafix 2016, C3
    
    core_name = []
    Lat = []
    Long = []
    for i in range(0,len(cd)):
        core_name += [cd['U.S. Geological Survey site name'][i]]
        Lat += [convert_coor_to_dec(cd['Latitude'][i])]
        Long += [convert_coor_to_dec(cd['Longitude'][i])]
        print (cd['U.S. Geological Survey site name'][i],cd['Latitude'][i],cd['Longitude'][i],\
            convert_coor_to_dec(cd['Latitude'][i]),convert_coor_to_dec(cd['Longitude'][i]))

    df = pd.DataFrame(
        {
            "Core Name": core_name,
            "Latitude": Lat,
            "Longitude": Long,
        }
    )

    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.Longitude, df.Latitude), crs="EPSG:4326")

    gdf.to_file(core_shpfile, driver='ESRI Shapefile')


import re
def convert_coor_to_dec(coor):
    deg, minutes, seconds, direction =  re.split('[°’”]', coor)
    return(float(deg) + float(minutes)/60 + float(seconds)/(60*60)) * (-1 if direction in [' W', ' S', ' W', ' S'] else 1)

if __name__ == "__main__":
    core_data = snakemake.input["core_data"]
    core_shpfile = snakemake.output["core_shpfile"]
    main(
        core_data,
        core_shpfile
    )