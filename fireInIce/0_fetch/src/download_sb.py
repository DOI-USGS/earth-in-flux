import os
import sciencebasepy

def main(sb_id,zipfile):
    # make out_dir if it doesn't exist
    out_dir = os.path.dirname(zipfile)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    # download the file
    sb = sciencebasepy.SbSession()
    sb_item = sb.get_item(sb_id)
    sb.get_item_files_zip(sb_item,out_dir)

if __name__ == "__main__":
    sb_id = snakemake.params["sb_id"]
    zipfile = snakemake.output["zipfile"]
    main(sb_id,zipfile)