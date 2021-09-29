import os

from snakemake.shell import shell

# Get rule information (input/output/params...)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

options = snakemake.params.get("options", "")

shell("""bamtools index -in {snakemake.input[0]} {options}  {log}""")

