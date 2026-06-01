CMD = """\
fastq_files=({input.fastq})
outdir=$(dirname {output.html})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    input_opt="--input-file1 ${{fastq_files[0]}} --input-file2 ${{fastq_files[1]}}"
    unc_name="unclassified#.fastq"
else
    input_opt="--input-file1 ${{fastq_files[0]}}"
    unc_name="unclassified.fastq"
fi
if [ "{params.store_unclassified}" = "True" ]; then
    unc_opt="--unclassified-out $unc_name"
else
    unc_opt=""
fi
sequana_taxonomy $input_opt \
    --databases {params.databases} \
    --output-directory $outdir \
    --thread {threads} \
    $unc_opt {params.options} > {log} 2>&1
if [ "{params.store_unclassified}" = "True" ] && [ ${{#fastq_files[@]}} -eq 2 ]; then
    if [ -f $outdir/kraken/unclassified_1.fastq ]; then
        r1=$outdir/kraken/unclassified_1.fastq
        r2=$outdir/kraken/unclassified_2.fastq
    else
        r1=$outdir/kraken/unclassified1.fastq
        r2=$outdir/kraken/unclassified2.fastq
    fi
    cat $r1 $r2 > $outdir/kraken/unclassified.fastq 2>> {log}
    rm -f $r1 $r2
elif [ "{params.store_unclassified}" != "True" ]; then
    touch $outdir/kraken/unclassified.fastq
fi
"""
