import re

# remove annoying metadata that sets vue warnings off
def remove_metadata(infile, outfile):
    output = open(outfile, "w")
    input = open(infile).read()
    output.write(re.sub("<metadata>.*?</metadata>\n", "", input, flags=re.DOTALL))
    output.close()

if __name__ == "__main__":
    infile = snakemake.input["infile"]
    outfile = snakemake.output["outfile"]
    remove_metadata(infile,outfile)