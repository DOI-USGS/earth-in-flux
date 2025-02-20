
#code modified from https://stackoverflow.com/questions/10454316/how-to-project-and-resample-a-grid-to-match-another-grid-with-gdal-python
from osgeo import gdal, gdalconst
import os

def main(src_file,match_file,dst_file):
    out_dir = os.path.dirname(dst_file)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Source
    src = gdal.Open(src_file, gdalconst.GA_ReadOnly)
    src_proj = src.GetProjection()

    # We want a section of source that matches this:
    match_ds = gdal.Open(match_file, gdalconst.GA_ReadOnly)
    match_proj = match_ds.GetProjection()
    match_geotrans = match_ds.GetGeoTransform()
    wide = match_ds.RasterXSize
    high = match_ds.RasterYSize

    # Output / destination
    dst = gdal.GetDriverByName('GTiff').Create(dst_file, wide, high, 1, gdalconst.GDT_Float32)
    dst.SetGeoTransform( match_geotrans )
    dst.SetProjection( match_proj)

    # Do the work
    gdal.ReprojectImage(src, dst, src_proj, match_proj, gdalconst.GRA_Bilinear)

    del dst # Flush

if __name__ == "__main__":
    src_file = snakemake.input["src_file"]
    match_file = snakemake.input["match_file"]
    dst_file = snakemake.output["dst_file"]
    main(src_file,match_file,dst_file)