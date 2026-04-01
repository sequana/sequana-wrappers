CMD = """\
bwa mem -t {threads} {params.options} {input.reference} {input.fastq} \
    > {output.sam} 2> {log}
sambamba view -t {threads} -S -f bam -o {output.bam} {output.sam} >> {log} 2>&1
sambamba sort {output.bam} -o {output.sorted} -t {threads} >> {log} 2>&1
"""
