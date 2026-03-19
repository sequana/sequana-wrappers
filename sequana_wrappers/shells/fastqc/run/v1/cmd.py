CMD = """\
mkdir -p {params.working_directory}
for fastq in {input}; do
    if [[ "$fastq" == *.bam ]] || [[ "$fastq" == *.sam ]]; then
        fastqc -t {threads} --outdir {params.working_directory} "$fastq" {params.options} >> {log} 2>&1
    else
        fastqc -t {threads} --outdir {params.working_directory} -f fastq "$fastq" {params.options} >> {log} 2>&1
    fi
done
touch {output.done}
"""
