CMD = """\
fastq_files=({input.fastq})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    reads="--in1 ${{fastq_files[0]}} --in2 ${{fastq_files[1]}}"
    out2_opt="--out2 {params.out2}"
else
    reads="--in1 ${{fastq_files[0]}}"
    out2_opt=""
fi
fastp --thread {threads} {params.options} {params.adapters} \
    $reads --out1 {output.r1} $out2_opt \
    --html {output.html} --json {output.json} > {log} 2>&1
"""
