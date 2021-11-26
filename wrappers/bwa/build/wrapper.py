__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Dev Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD"


from snakemake.shell import shell


reference = snakemake.input.get("reference", "")
options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
index = snakemake.params.get("index_algorithm")
assert index in ['is', 'bwtsw', 'rb2'], "index must be set to 'is', 'rb2' or bwtsw'; please see bwa documentation"


shell("bwa index -a {index} {options} {reference}  {log} ")

shell("samtools faidx {reference}")
