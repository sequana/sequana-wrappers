CMD = """\
sambamba view {params.options} \
    --format=bam \
    --filter='mapping_quality >= {params.threshold}' \
    -o {output.bam} \
    {input.bam} > {log} 2>&1
"""
