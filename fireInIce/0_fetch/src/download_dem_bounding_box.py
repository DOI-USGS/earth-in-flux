import os
import numpy as np
import geopandas as gpd
from bmi_topography import Topography
from shapely.geometry import Polygon

def snakemake_type_exists(snakemake_type,string,default_input):
    if hasattr(snakemake_type, string):
        return snakemake_type[string]
    else:
        return default_input

def bmi_download_dem(UL_corner, LR_corner,extent_shpfile, demfile, data_product, dem_product, buffer):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(demfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # set parameters of dem download
    params = Topography.DEFAULT.copy()
    params["data_type"] = data_product
    params["dem_type"] = dem_product
    params["output_format"] = "GTiff"
    params["west"] = float(UL_corner[1])
    params["south"] = float(LR_corner[0])
    params["east"] = float(LR_corner[1])
    params["north"] = float(UL_corner[0])
    params["cache_dir"] = out_dir

    # set api key if exists (you need this when you run out of free requests)
    if os.path.isfile("opentopography_api_key.txt"):
        with open("opentopography_api_key.txt", "r") as f:
            api_key = f.read()
        params["api_key"] = api_key

    # setup request and fetch
    topo_request = Topography(**params)
    topo_request.fetch()

    # rename the file 
    os.rename(
        out_dir
        + "/"
        + str(params["dem_type"])
        + "_"
        + str(params["south"])
        + "_"
        + str(params["west"])
        + "_"
        + str(params["north"])
        + "_"
        + str(params["east"])
        + ".tif",
        demfile,
    )

    domain_geom = Polygon(
        zip(
            [UL_corner[1],UL_corner[1], LR_corner[1],LR_corner[1]],
            [UL_corner[0], LR_corner[0], LR_corner[0], UL_corner[0]],
        )
    )
    domain_polygon = gpd.GeoDataFrame(index=[0], crs="EPSG:4326", geometry=[domain_geom])
    domain_polygon.to_file(extent_shpfile)


if __name__ == "__main__":
    data_product = snakemake_type_exists(snakemake.params,"data_product","/API/globaldem")
    dem_product = snakemake_type_exists(snakemake.params,"dem_product","NASADEM")
    buffer = snakemake.params["buffer"]
    UL_corner = snakemake.params["UL_corner"]
    LR_corner = snakemake.params["LR_corner"]
    extent_shpfile = snakemake.output["extent_shpfile"]
    demfile = snakemake.output["demfile"]
    download_mode = snakemake_type_exists(snakemake.params,"download_mode","bmi")
    bmi_download_dem(UL_corner, LR_corner, extent_shpfile, demfile, data_product, dem_product, buffer)
