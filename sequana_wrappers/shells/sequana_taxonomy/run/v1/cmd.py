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
if [ "{params.store_unclassified}" = "True" ]; then
    if [ "${{#fastq_files[@]}}" -eq 2 ]; then
        cmd="$cmd --unclassified-out unclassified#.fastq"
    else
        cmd="$cmd --unclassified-out unclassified.fastq"
    fi
fi
eval "$cmd" > {log} 2>&1
if [ "{params.store_unclassified}" = "True" ]; then
    unclass_dir="$outdir/kraken"
    if [ -f "$unclass_dir/unclassified_1.fastq" ] && [ -f "$unclass_dir/unclassified_2.fastq" ]; then
        cat "$unclass_dir/unclassified_1.fastq" "$unclass_dir/unclassified_2.fastq" > "$unclass_dir/unclassified.fastq"
    elif [ ! -f "$unclass_dir/unclassified.fastq" ]; then
        touch "$unclass_dir/unclassified.fastq"
    fi
fi
"""
