__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Dev Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD"


from snakemake.shell import shell


fastq = snakemake.input.fastq
reference = snakemake.input.reference
options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
output_sorted_bam = snakemake.output.sorted
params = snakemake.params


from easydev import cmd_exists
if cmd_exists("pbwa"):
    bwa_exe = "pbwa"
else:
    bwa_exe = "bwa"


sambamba_sort = ""

shell("""
        ({bwa_exe} mem -t {snakemake.threads} {options} \
        {reference} {fastq} | \
        sambamba view -t {snakemake.threads} -S -f bam -o /dev/stdout /dev/stdin | \
        sambamba sort /dev/stdin -o {output_sorted_bam} -t {snakemake.threads} \
        --tmpdir={params.tmp_directory} """ + sambamba_sort + """ )  {log}
"""
)


