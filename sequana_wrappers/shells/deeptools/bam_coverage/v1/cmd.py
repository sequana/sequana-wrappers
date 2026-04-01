CMD = """\
if [ ! -f {input.bam}.bai ]; then
    samtools index {input.bam} >> {log} 2>&1
fi
bamCoverage --bam {input.bam} \
    --outFileName {output} \
    --numberOfProcessors {threads} \
    {params.options} > {log} 2>&1
"""
