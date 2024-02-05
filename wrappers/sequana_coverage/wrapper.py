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

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get rule information (input/output/params...)
# inputs
input_bed = snakemake.input.bed
input_fasta = snakemake.input.fasta

# outputs
output_html = snakemake.output[0]
report_dir = output_html.rsplit("/", 1)[0]

# default values from sequana_coverage
annotation = snakemake.params.get("annotation", None)
circular = snakemake.params.get("circular", True)
chunksize = snakemake.params.get("chunksize", 5000000)
double_threshold = snakemake.params.get("double_threshold", 0.5)
gc_window_size = snakemake.params.get("gc_window_size", 101)
high = snakemake.params.get("high_threshold", 4)
low = snakemake.params.get("low_threshold", -4)
params_k = snakemake.params.get("mixture_models", 2)
options = snakemake.params.get("options", "")
window_size = snakemake.params.get("window_size", 20001)
output_directory = snakemake.params.get("output_directory", "report")

# create the command
cmd = f"sequana_coverage --input-file {input_bed} -H {high} -L {low} "\
       "--clustering-parameter {double_threshold} --chunk-size {chunksize} "\
       "--window-gc {gc_window_size} --mixture-models {params_k} "\
       "--output-directory {output_directory} --window-median {window_size}"

if circular:
    cmd += " -o "

if annotation:
    cmd += f" --annotation-file {annotation} "

if input_fasta:
    cmd += f" --reference-file {input_fasta} "

cmd += " {options} {log} "

shell(cmd)
