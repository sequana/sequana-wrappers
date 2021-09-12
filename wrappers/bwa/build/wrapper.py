__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Dev Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD"


from snakemake.shell import shell


reference = snakemake.input.get("reference", "")
options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
index = snakemake.options.index_algorithm


from easydev import cmd_exists
if cmd_exists("pbwa"):
    cmd = "pbwa"
else:
    cmd = "bwa"


shell(
    "{cmd} index -a {index} {options} {reference} > {log} "
    "samtools faidx {reference}"
)
