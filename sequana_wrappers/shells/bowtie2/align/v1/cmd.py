CMD = """\
(
  fastq_files=({input.fastq})
  if [ ${{#fastq_files[@]}} -eq 2 ]; then
    reads="-1 ${{fastq_files[0]}} -2 ${{fastq_files[1]}}"
  else
    reads="-U ${{fastq_files[0]}}"
  fi
  idx="{input.idx[0]}"
  idx="${{idx%.1.bt2}}"
  idx="${{idx%.1.bt2l}}"
  bowtie2 --threads {threads} {params.options} -x $idx $reads \
  | samtools sort - > {output.bam} \
  && samtools index {output.bam}
) > {log} 2>&1
"""
