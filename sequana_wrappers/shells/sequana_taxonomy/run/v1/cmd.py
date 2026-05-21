CMD = """\
fastq_files=({input.fastq})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    input_opt="--input-file1 ${{fastq_files[0]}} --input-file2 ${{fastq_files[1]}}"
else
    input_opt="--input-file1 ${{fastq_files[0]}}"
fi
outdir=$(dirname {output.html})
cmd="sequana_taxonomy $input_opt \
    --databases {params.databases} \
    --output-directory $outdir \
    --thread {threads} \
    {params.options}"
eval "$cmd" > {log} 2>&1
if [ "{params.store_unclassified}" = "True" ]; then
    touch $outdir/unclassified.fastq
fi
"""
