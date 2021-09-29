import os

from snakemake.shell import shell

# Get rule information (input/output/params...)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("""bamtools sort -in {snakemake.input[0]} -out {snakemake.output[0]} {snakemake.params.options} {log}""")

