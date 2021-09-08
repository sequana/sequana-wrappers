__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Dev Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD"


from snakemake.shell import shell


prefix = snakemake.input.index.replace(".1.ebwt", "")

options = snakemake.params.get("options", "")
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# fastq could be a single filename or a list
reads = snakemake.input.fastq

print('--------')
print(snakemake.input.fastq)
print(reads)

try:
    # if a single string is provided, this is a Single-End filename
    # we just use strip to fall back 
    reads = reads.strip()
    reads = f" {reads}"
except AttributeError:
    reads = f"-1 {reads[0]} -2 {reads[1]}"


print(reads)

shell(
    "(bowtie -S {options} -p {threads} -x {prefix} {reads}"
    "| samtools view -Sbh  - > {snakemake.output.bam}) {log}"
)


# sort result
shell("samtools sort -o {snakemake.output.sort} {snakemake.output.bam}")
shell("samtools index {snakemake.output.sort}")

