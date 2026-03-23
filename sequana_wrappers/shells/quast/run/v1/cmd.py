CMD = """\
fastq_files=({input.fastq})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    input_fastq="-1 ${{fastq_files[0]}} -2 ${{fastq_files[1]}}"
else
    input_fastq="--single ${{fastq_files[0]}}"
fi
outdir=$(dirname {output})
quast.py {params.options} \
    -t {threads} \
    {input.assembly} \
    $input_fastq \
    --output-dir $outdir > {log} 2>&1
touch {output}
"""
