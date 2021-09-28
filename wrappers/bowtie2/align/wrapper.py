
from snakemake.shell import shell

options = snakemake.params.get("options", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

n = len(snakemake.input.fastq)

if n == 1:
    reads = "-U {}".format(*snakemake.input.fastq)
elif n == 2:
    reads = "-1 {} -2 {}".format(*snakemake.input.fastq)
else:
    raise ValueError("input->fastq must have 1 (single-end) or 2 (paired-end) elements.")

shell(
    "(bowtie2 --threads {snakemake.threads} {options} "
    "-x {snakemake.params.index} {reads} "
    "| samtools view -Sbh -o {snakemake.output.bam} -) {log}"
)

# index the bam file    
shell("bamtools index -in {snakemake.output.bam}")

# sort the bam
shell("bamtools sort -in {snakemake.output.bam} -out {snakemake.output.sorted}")

# and index it
shell("bamtools index -in {snakemake.output.sorted}")
