import os
import zipfile as zf
import py7zr

def unzip_file(zipfile, out_dir):
    if zipfile.endswith(".7z"):
        with py7zr.SevenZipFile(zipfile) as zip_ref:
            zip_ref.extractall(out_dir)
    elif zipfile.endswith (".zip"):
        with zf.ZipFile(zipfile, "r") as zip_ref:
            zip_ref.extractall(out_dir)

def main(zipfile):
    # make out_dir if it doesn't exist
    tmp_dir = os.path.dirname(zipfile)
    if zipfile.endswith(".7z"):
        out_dir = tmp_dir + '/' + os.path.basename(zipfile[:-3])
    elif zipfile.endswith (".zip"):
        out_dir = tmp_dir + '/' + os.path.basename(zipfile[:-4])
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # unzip the file
    unzip_file(zipfile, out_dir)


if __name__ == "__main__":
    zipfile = snakemake.input["zipfile"]
    main(zipfile)
