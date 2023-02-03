"""Snakemake wrapper for digital normalisation"""
from pathlib import Path

from snakemake.shell import shell


# we want to append the log
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


options = snakemake.params.get("options", "")
params = snakemake.params
threads = snakemake.threads
output = snakemake.output

graph = snakemake.params.get("tmp_graph_name", "graph.ct")


# Files name without .gz extension
fastq = snakemake.input.get("fastq")
if isinstance(fastq, str):
    paired = False
    # use list even for single-end data
    fastq = [fastq]
elif len(fastq) == 2:
    paired = True
else:
    raise ValueError("digital_normalisation wrapper expects a list of 2 fastq files for paired data, or a string for single-end data")
    sys.exit(1)

# Check if FASTQ are compressed. If so,
fastq_to_delete = []

for filename in fastq:
    if filename.endswith(".gz"):
        shell(f"unpigz -p {threads} -fk {filename}")
        fastq_to_delete.append(filename.strip(".gz"))


if paired:
    prefix = Path(output.fastq[0]).stem

    shell("interleave-reads.py {fastq[0]} {fastq[1]} --output {prefix}.pe  {log} ")

    # Digital normalisation
    shell("""normalize-by-median.py --paired --ksize {params.ksize} \
    --cutoff {params.cutoff} -M {params.m} {params.options} \
    --savegraph {graph} {prefix}.pe \
    --output {prefix}.pe.keep  {log}
    """)

    # Filter abundance
    shell("""filter-abund.py --threads {threads} -V {graph} \
        {prefix}.pe.keep --output {prefix}.pe.filter \
         {log} """)

    # Extract paired reads
    shell("""extract-paired-reads.py {prefix}.pe.filter \
        --output-paired {prefix}.pe \
        --output-single {prefix}.orphans  {log}
    """)

    # Split paired reads
    shell("""split-paired-reads.py {prefix}.pe \
        -1 {output.fastq[0]} -2 {output.fastq[1]}  {log}
    """)
else:
    prefix = Path(output.fastq).stem
    # Digital normalisation
    shell("""normalize-by-median.py  --ksize {params.ksize} \
        --cutoff {params.cutoff} -M {params.m} {params.options} \
        --savegraph {graph} {fastq[0]} \
        --output {prefix}.se.keep  {log}
    """)

    # Filter abundance
    shell("""filter-abund.py --threads {threads} -V {graph} \
        {prefix}.se.keep --output {prefix}.se.filter \
         {log} """)

    shell("""mv {prefix}.se.filter {output.fastq}""")


# cleanup
for filename in fastq_to_delete:
    shell(f"rm -f {filename}")
shell(f"rm -f {graph}")
if paired:
    shell(f"rm -f {prefix}.pe {prefix}.pe.filter {prefix}.pe.keep {prefix}.orphans")
else:
    shell(f"rm -f {prefix}.se.keep {prefix}.se.filter")
