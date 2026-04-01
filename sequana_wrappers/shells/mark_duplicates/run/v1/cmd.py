CMD = """\
mkdir -p {params.tmpdir}
picard MarkDuplicates \
    I={input.bam} \
    O={output.bam} \
    M={output.metrics} \
    REMOVE_DUPLICATES={params.remove_dup} \
    TMP_DIR={params.tmpdir} \
    {params.options} > {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
"""
