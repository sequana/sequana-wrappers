__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import os

from snakemake.shell import shell


# Get rule information (input/output/params...)
input_samplesheet = snakemake.input['samplesheet']

# figure out the expected output directory
output_json = snakemake.output[0]

# params and options
params = snakemake.params
options = snakemake.params.get("options", "")


################### The code 


cmd = "bcl2fastq -p {snakemake.threads} --barcode-mismatches {params.barcode_mismatch}"
cmd += " --runfolder-dir {params.indir}"
cmd += " --intensities-dir {params.indir}/Data/Intensities"
if input_samplesheet.strip()!= "":
    cmd += " --sample-sheet {input_samplesheet}"
cmd += " --output-dir {}".format(os.path.abspath("."))

if params.ignore_missing_bcls:
    cmd += " --ignore-missing-bcls "

if params.no_bgzf_compression:
    cmd += " --no-bgzf-compression "

if params.merge_all_lanes:
    cmd += " --no-lane-splitting "

cmd += options

shell(cmd)
