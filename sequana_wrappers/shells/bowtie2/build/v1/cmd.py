CMD = """\
idx="{output[0]}"
idx="${{idx%.1.bt2}}"
idx="${{idx%.1.bt2l}}"
bowtie2-build --threads {threads} {params.options} {input.reference} $idx > {log} 2>&1
"""
