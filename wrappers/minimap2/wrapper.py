__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import os

from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fastq = snakemake.input['fastq']
input_reference = snakemake.input['reference']
output = snakemake.output
params = snakemake.params


shell("""
minimap2 -t {threads} {input_reference} {input_fastq} {params.options} -a | samtools view -b | bamtools sort -in - -out {output}
""")
