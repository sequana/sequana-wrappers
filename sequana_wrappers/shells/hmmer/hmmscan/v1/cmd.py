CMD = """\
profile="{input.profile}"
profile="${{profile%.h3f}}"
profile="${{profile%.h3i}}"
profile="${{profile%.h3m}}"
profile="${{profile%.h3p}}"
hmmscan {params.options} \
    --cpu {threads} \
    --tblout {output.tblout} \
    --domtblout {output.domtblout} \
    $profile {input.fasta} > {log} 2>&1
"""
