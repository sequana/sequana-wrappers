CMD = """\
idx="{output[0]}"
idx="${{idx%.1.ebwt}}"
bowtie-build --threads {threads} {params.options} {input.reference} $idx > {log} 2>&1
samtools faidx {input.reference} >> {log} 2>&1
"""
