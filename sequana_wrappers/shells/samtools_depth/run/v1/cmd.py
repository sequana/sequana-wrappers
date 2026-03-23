CMD = """\
samtools depth -m 20000 -aa {input.bam} {params.options} > {output} 2>{log}
"""
