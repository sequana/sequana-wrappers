CMD = """\
output_file="{output[0]}"
output_dir=$(dirname "$output_file")
prefix=$(basename "$output_file" .contigs.fasta)
rm -f $output_dir/canu.done $output_dir/canu.failed
canu {params.step} \
    -p $prefix \
    -d $output_dir \
    genomeSize={params.genome_size} \
    maxThreads={threads} \
    {params.options} \
    -{params.preset} {input} > {log} 2>&1
"""
