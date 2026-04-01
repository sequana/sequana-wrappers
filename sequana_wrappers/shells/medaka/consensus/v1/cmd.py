CMD = """\
output_dir=$(dirname {output.fasta})
default_file=$output_dir/consensus.fasta
(medaka_consensus {params.options} \
    -t {threads} \
    -i {input.fastq} \
    -d {input.assembly} \
    -o $output_dir \
    -m {params.model} \
    && mv $default_file {output.fasta}) > {log} 2>&1
"""
