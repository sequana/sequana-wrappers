__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Dev Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD"


from snakemake.shell import shell

# The input reference in fasta format
reference = snakemake.input['reference']

# prefix must be extracted
prefix = snakemake.output[0].replace(".1.ebwt", "")

# options is optional
options = snakemake.params.get("options", "")

# the threading and logging
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


shell("bowtie-build {reference} {prefix} --threads {threads} {options} {log}")
shell("samtools faidx {reference}  {log}")
