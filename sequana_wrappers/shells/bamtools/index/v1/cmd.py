CMD = """\
bamtools index -in {input.bam} {params.options} > {log} 2>&1
"""
