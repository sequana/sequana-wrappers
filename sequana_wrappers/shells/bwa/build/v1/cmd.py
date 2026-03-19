CMD = """\
bwa index -a {params.index_algorithm} {params.options} {input.reference} > {log} 2>&1
samtools faidx {input.reference} 2>> {log}
"""
