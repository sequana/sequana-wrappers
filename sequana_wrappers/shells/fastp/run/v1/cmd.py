CMD = """\
fastq_files=({input.fastq})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    reads="--in1 ${{fastq_files[0]}} --in2 ${{fastq_files[1]}}"
    trimmed_opt="--out1 {output.r1} --out2 {output.r2}"
else
    reads="--in1 ${{fastq_files[0]}}"
    trimmed_opt="--out1 {output.r1}"
fi
fastp --thread {threads} {params.options} {params.adapters} \
    $reads $trimmed_opt \
    --html {output.html} --json {output.json} > {log} 2>&1
"""
