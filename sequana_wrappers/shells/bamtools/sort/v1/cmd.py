CMD = """\
bamtools sort -in {input.bam} -out {output.bam} {params.options} > {log} 2>&1
"""
