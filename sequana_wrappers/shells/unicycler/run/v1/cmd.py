CMD = """\
fastq_files=({input.fastq})
if [ "{params.long_reads}" = "True" ]; then
    input_opt="-l ${{fastq_files[0]}}"
elif [ ${{#fastq_files[@]}} -eq 2 ]; then
    input_opt="-1 ${{fastq_files[0]}} -2 ${{fastq_files[1]}}"
else
    input_opt="-s ${{fastq_files[0]}}"
fi
outdir=$(dirname {output.fasta})
unicycler --mode {params.mode} \
    --threads {threads} \
    $input_opt \
    -o $outdir \
    {params.options} > {log} 2>&1
cp $outdir/assembly.fasta {output.fasta}
"""
