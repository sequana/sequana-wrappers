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
    awk '$1=="U"{{print $2}}' $outdir/kraken/kraken.out > $outdir/kraken/unclassified_ids.txt
    if [ ${{#fastq_files[@]}} -eq 2 ]; then
        seqkit grep -f $outdir/kraken/unclassified_ids.txt ${{fastq_files[0]}} > $outdir/kraken/_unclassified_R1.fastq 2>> {log}
        seqkit grep -f $outdir/kraken/unclassified_ids.txt ${{fastq_files[1]}} > $outdir/kraken/_unclassified_R2.fastq 2>> {log}
        cat $outdir/kraken/_unclassified_R1.fastq $outdir/kraken/_unclassified_R2.fastq > $outdir/kraken/unclassified.fastq
        rm -f $outdir/kraken/_unclassified_R1.fastq $outdir/kraken/_unclassified_R2.fastq
    else
        seqkit grep -f $outdir/kraken/unclassified_ids.txt ${{fastq_files[0]}} > $outdir/kraken/unclassified.fastq 2>> {log}
    fi
    rm -f $outdir/kraken/unclassified_ids.txt
else
    touch $outdir/kraken/unclassified.fastq
fi
"""
