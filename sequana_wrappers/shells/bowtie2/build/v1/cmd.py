CMD = """\
bowtie2-build --threads {threads} {params.options} {input.reference} reference/genome > {log} 2>&1
"""
