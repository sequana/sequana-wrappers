#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Documentation: http://sequana.readthedocs.io
#  Contributors:  https://github.com/sequana/sequana/graphs/contributors
##############################################################################
import os

from snakemake.shell import shell

# Get rule information (input/output/params...)
input_samplesheet = snakemake.input["samplesheet"]

# figure out the expected output directory
output_json = snakemake.output[0]

# params and options
params = snakemake.params
options = snakemake.params.get("options", "")


################### The code

cmd = "bcl2fastq -p {snakemake.threads} --barcode-mismatches {params.barcode_mismatch}"
cmd += " --runfolder-dir {params.indir}"
cmd += " --intensities-dir {params.indir}/Data/Intensities"
if input_samplesheet.strip() != "":
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
