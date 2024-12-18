
import os
import numpy as np
from osgeo import gdal, osr

def calc_bedrock(dem_file,ice_thickness_file,ice_height_file,bedrock_file):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(bedrock_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # read in the dem
    ds = gdal.Open(dem_file)
    ds_it = gdal.Open(ice_thickness_file)

    #get metadata
    gt = ds.GetGeoTransform()
    cs = ds.GetProjection()
    dtype = ds.GetRasterBand(1).DataType
    cellsx = ds.GetRasterBand(1).XSize
    cellsy = ds.GetRasterBand(1).YSize
    NaN = ds.GetRasterBand(1).GetNoDataValue()
    if NaN == None:
        NaN = -9999

    #convert to numpy array
    dem = np.array(ds.GetRasterBand(1).ReadAsArray())
    it = np.array(ds_it.GetRasterBand(1).ReadAsArray())

    #calc bedrock
    bedrock = dem - it

    #ice_height
    ih = np.zeros_like(dem) - 100000.
    ih[it!=0] = dem[it!=0]
    
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(ice_height_file,cellsx,cellsy, 1, gdal.GDT_Float32)
    dataset.GetRasterBand(1).WriteArray(ih)
    dataset.GetRasterBand(1).SetNoDataValue(NaN)

    dataset.SetGeoTransform(gt)
    dataset.SetProjection(cs)
    dataset.FlushCache()
    dataset=None

    dataset = driver.Create(bedrock_file,cellsx,cellsy, 1, gdal.GDT_Float32)
    dataset.GetRasterBand(1).WriteArray(bedrock)
    dataset.GetRasterBand(1).SetNoDataValue(NaN)

    dataset.SetGeoTransform(gt)
    dataset.SetProjection(cs)
    dataset.FlushCache()
    dataset=None

if __name__ == "__main__":
    dem_file = snakemake.input["dem_file"]
    ice_thickness_file = snakemake.input["ice_thickness_file"]
    ice_height_file = snakemake.output["ice_height_file"]
    bedrock_file = snakemake.output["bedrock_file"]
    calc_bedrock(dem_file,ice_thickness_file,ice_height_file,bedrock_file)