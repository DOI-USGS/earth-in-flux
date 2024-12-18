from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
import os
import numpy as np
from osgeo import gdal, osr

x_nei = [-1,0,1,-1,1,-1,0,1]
y_nei = [1,1,1,0,0,-1,-1,-1]

def seed_ice_hole(it,cellsx,cellsy,it_threshold):
    for x in range(0,cellsx):
        for y in range(0,cellsy):
            if it[y,x] > it_threshold:
                for nei in range(0,8):
                    if it[y+y_nei[nei],x+x_nei[nei]] == 0:
                        it[y+y_nei[nei],x+x_nei[nei]] = -9999
    return it
            
def fill_ice_hole(it,cellsx,cellsy):
    trigger = 0
    fill_iteration = 0
    while trigger == 0:
        trigger = 1
        for x in range(0,cellsx):
            for y in range(0,cellsy):
                if it[y,x] == -9999:
                    trigger = 0
                    mini_trigger = 0
                    for nei in range(0,8):
                        if it[y+y_nei[nei],x+x_nei[nei]] == 0:
                            mini_trigger = 1
                            it[y+y_nei[nei],x+x_nei[nei]] = -9999
                    if mini_trigger == 0:
                        it[y,x] = -999999
        fill_iteration += 1
        print ("filling ice, try " + str(fill_iteration) + "...")
    return it

def remove_nans(ice_thickness_file,fixed_ice_thickness_file,convolution_window):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(fixed_ice_thickness_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # read in the dem
    ds_it = gdal.Open(ice_thickness_file)

    #get metadata
    gt = ds_it.GetGeoTransform()
    cs = ds_it.GetProjection()
    dtype = ds_it.GetRasterBand(1).DataType
    cellsx = ds_it.GetRasterBand(1).XSize
    cellsy = ds_it.GetRasterBand(1).YSize
    NaN =0

    #convert to numpy array
    it = np.array(ds_it.GetRasterBand(1).ReadAsArray())
    
    #convert zero to NaN
    it = seed_ice_hole(it,cellsx,cellsy,750.)
    it = fill_ice_hole(it,cellsx,cellsy)
    it[it==-999999] = np.nan

    #set kernal to one neighbor
    kernel = Gaussian2DKernel(x_stddev=convolution_window)

    # create a "fixed" image with NaNs replaced by interpolated values
    fixed_it= interpolate_replace_nans(it, kernel)

    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(fixed_ice_thickness_file,cellsx,cellsy, 1, gdal.GDT_Float32)
    dataset.GetRasterBand(1).WriteArray(fixed_it)
    dataset.GetRasterBand(1).SetNoDataValue(NaN)
    dataset.SetGeoTransform(gt)
    dataset.SetProjection(cs)
    dataset.FlushCache()
    dataset=None

if __name__ == "__main__":
    convolution_window = snakemake.params["convolution_window"]
    ice_thickness_file = snakemake.input["ice_thickness_file"]
    fixed_ice_thickness_file = snakemake.output["ice_thickness_file"]
    remove_nans(ice_thickness_file,fixed_ice_thickness_file,convolution_window)