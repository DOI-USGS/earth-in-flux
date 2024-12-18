import os
import pandas as pd
import geopandas as gpd
from pygeohydro import WBD

def main(huc_list, extent_shpfile):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(extent_shpfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #download huc
    wbd_list = []
    for huc in huc_list:
        huc_level = str(len(huc))
        wbd = WBD("huc" + huc_level)
        wbd_list += [wbd.byids("huc" + huc_level, [huc]).dissolve()]
    
    extent = gpd.GeoDataFrame(pd.concat(wbd_list)).dissolve()
    extent.to_file(extent_shpfile)

if __name__ == "__main__":
    huc_list = snakemake.params["huc_list"]
    extent_shpfile = snakemake.output["extent_shpfile"]
    main(huc_list, extent_shpfile)
