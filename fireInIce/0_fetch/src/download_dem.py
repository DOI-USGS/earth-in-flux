import os
import geopandas as gpd
from bmi_topography import Topography
import py3dep
from shapely.geometry import Polygon

def snakemake_type_exists(snakemake_type,string,default_input):
    if hasattr(snakemake_type, string):
        return snakemake_type[string]
    else:
        return default_input

def get_shapefile_extent(shapefile):
    west, south, east, north = shapefile.geometry.total_bounds
    length = east - west
    height = north - south
    return west, east, length, south, north, height


def buffer_shapefile(west, east, length, south, north, height, buffer):
    west -= buffer / 100.0 * length
    south -= buffer / 100.0 * height
    east += buffer / 100.0 * length
    north += buffer / 100.0 * height
    return west, east, south, north


def bmi_download_dem(extent_shpfile, demfile, data_product, dem_product, buffer):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(demfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load shapefile
    shp = gpd.read_file(extent_shpfile)
    # get extent of shapefile #assums WGS84 or some sort of geographic coordinate system, Potential FIX: could convert if necessary.
    west, east, length, south, north, height = get_shapefile_extent(shp)
    buff_west, buff_east, buff_south, buff_north = buffer_shapefile(
        west, east, length, south, north, height, buffer
    )

    # set parameters of dem download
    params = Topography.DEFAULT.copy()
    params["data_type"] = data_product
    params["dem_type"] = dem_product
    params["output_format"] = "GTiff"
    params["west"] = buff_west
    params["south"] = buff_south
    params["east"] = buff_east
    params["north"] = buff_north
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
        out_dir + "/dem.tif",
    )

def py3dep_download_dem(extent_shpfile, demfile, py3dep_resolution, buffer):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(demfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # load shapefile
    shp = gpd.read_file(extent_shpfile)
    # get extent of shapefile #assums WGS84 or some sort of geographic coordinate system, Potential FIX: could convert if necessary.
    west, east, length, south, north, height = get_shapefile_extent(shp)
    buff_west, buff_east, buff_south, buff_north = buffer_shapefile(
        west, east, length, south, north, height, buffer
    )

    #specify geometry to download dem from
    buffered_domain_polygon = Polygon(
        zip(
            [buff_west, buff_west, buff_east, buff_east],
            [buff_north, buff_south, buff_south, buff_north],
        )
    )

    #download dem
    dem = py3dep.get_map("DEM", buffered_domain_polygon, resolution = py3dep_resolution, geo_crs = shp.crs, crs = 4326)

    #save dem
    dem.rio.to_raster(out_dir + "/dem.tif")

if __name__ == "__main__":
    data_product = snakemake_type_exists(snakemake.params,"data_product","/API/globaldem")
    dem_product = snakemake_type_exists(snakemake.params,"dem_product","NASADEM")
    py3dep_resolution = snakemake_type_exists(snakemake.params,"py3dep_resolution",30)
    buffer = snakemake.params["buffer"]
    extent_shpfile = snakemake.input["extent_shpfile"]
    demfile = snakemake.output["demfile"]
    download_mode = snakemake_type_exists(snakemake.params,"download_mode","bmi")
    if download_mode == "bmi":
        bmi_download_dem(extent_shpfile, demfile, data_product, dem_product, buffer)
    elif download_mode == "py3dep":
        py3dep_download_dem(extent_shpfile, demfile, py3dep_resolution, buffer)
    else:
        raise Exception("Incorrect download mode. Please use bmi or py3dep.")
