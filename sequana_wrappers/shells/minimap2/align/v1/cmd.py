CMD = """\
minimap2 -t {threads} {input.reference} {input.fastq} {params.options} -a \
| samtools sort -@ {threads} -o {output}
"""
