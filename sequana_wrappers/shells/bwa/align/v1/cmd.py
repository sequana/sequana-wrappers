CMD = """\
mkdir -p {params.tmp_directory}
(bwa mem -t {threads} {params.options} {input.reference} {input.fastq} \
 | sambamba view -t {threads} -S -f bam -o /dev/stdout /dev/stdin \
 | sambamba sort /dev/stdin -o {output.sorted} -t {threads} --tmpdir={params.tmp_directory}) \
> {log} 2>&1
"""
