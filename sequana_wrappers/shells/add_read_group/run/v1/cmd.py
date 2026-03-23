CMD = """\
picard AddOrReplaceReadGroups \
    -VALIDATION_STRINGENCY SILENT \
    -I {input.bam} \
    -O {output.bam} \
    -PL {params.PL} \
    -LB {params.LB} \
    -PU {params.PU} \
    -SM {params.SM} \
    -ID {params.ID} \
    {params.options} > {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
"""
