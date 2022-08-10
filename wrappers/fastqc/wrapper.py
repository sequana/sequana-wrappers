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
input_fastq = snakemake.input

# figure out the expected output directory
done = snakemake.output[0]

# figure out the output directory
outdir = snakemake.params["working_directory"]

params = snakemake.params
# we do not use this alias because we will concatenate possibly 2 output of
# fastqc into a single log file for the paired data
# log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log = snakemake.log[0]


# Note that if the input file is empty, fastqc creates a HTML file
try:
    input_fastq = input_fastq.split()
except AttributeError:
    pass

os.makedirs(outdir, exist_ok=True)

for fastq_file in input_fastq:
    if fastq_file.endswith((".bam", "sam")):
        shell(" fastqc -t {snakemake.threads} --outdir {outdir} {fastq_file} {params.options} &>> {log}")
    else:
        shell(" fastqc -t {snakemake.threads} --outdir {outdir} -f fastq {fastq_file} {params.options} &>> {log}")

shell("touch {done} ")
