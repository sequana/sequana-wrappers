import os

from snakemake.shell import shell

from sequana.bedtools import GenomeCov
from sequana.modules_report.coverage import CoverageModule
from sequana.utils import config

# Get rule information (input/output/params...)


input_bam = snakemake.input.bam
input_ref = snakemake.input.ref

output_vcf = snakemake.output[0]

ploidy = snakemake.params['ploidy']
options = snakemake.params['options']

log = snakemake.log[0]


cmd = "freebayes {options} --ploidy {ploidy} -f {input_ref} -v {output_vcf} {input_bam} &>{log}"

shell(cmd)

