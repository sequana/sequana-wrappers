__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import os

from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fastq = snakemake.input

# figure out the expected output directory
done = snakemake.output[0]

# figure out the output directory 
outdir = snakemake.params['working_directory']

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

for fastq_file in input_fastq:
    if fastq_file.endswith((".bam", "sam")):
        shell(
        " fastqc -t {snakemake.threads} --outdir {outdir} "
            " {fast_file} {params.options} &>> {log}")
    else:
        shell(
            " fastqc -t {snakemake.threads} --outdir {outdir} -f fastq "
            " {fastq_file} {params.options} &>> {log}")

shell("touch {done} ")

