CMD = """\
(
  fastq_files=({input.fastq})
  if [ ${{#fastq_files[@]}} -eq 2 ]; then
    reads="-1 ${{fastq_files[0]}} -2 ${{fastq_files[1]}}"
  else
    reads="${{fastq_files[0]}}"
  fi
  idx="{input.index}"
  idx="${{idx%.1.ebwt}}"
  bowtie -S {params.options} -p {threads} -x $idx $reads \
  | samtools view -Sbh - > {output.bam}) > {log} 2>&1
samtools sort -@ {threads} -o {output.sorted} {output.bam} >> {log} 2>&1
samtools index {output.sorted} >> {log} 2>&1
"""
