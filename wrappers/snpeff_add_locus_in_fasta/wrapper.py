import os

from snakemake.shell import shell

# Get rule information (input/output/params...)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell("""bamtools sort -in {snakemake.input[0]} -out {snakemake.output[0]} {snakemake.params.options} {log}""")


from sequana import SnpEff
if input['ann'].endswith(".gbk"):
    snpeff = SnpEff(input['ann'], log=log['log'])
elif input['ann'].endswith("gff") or input['ann'].endswith('gff3'):
    snpeff = SnpEff(input['ann'], log=log['log'], fastafile=input['fasta'])
else:
    raise IOError("Your annotation file does not end with gbk or gff/gff3 extension")

snpeff.add_locus_in_fasta(input['fasta'], output['fasta'])
