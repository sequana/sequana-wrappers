CMD = """\
indexdir=$(dirname {input.index})
outdir=$(dirname {output.quant})
nfiles=$(echo {input.fastq} | wc -w)
if [ $nfiles -eq 2 ]; then
    R1=$(echo {input.fastq} | cut -d' ' -f1)
    R2=$(echo {input.fastq} | cut -d' ' -f2)
    salmon quant -i $indexdir {params.options} \
        -1 $R1 -2 $R2 -p {threads} \
        --validateMappings -o $outdir > {log} 2>&1
else
    salmon quant -i $indexdir {params.options} \
        -r {input.fastq} -p {threads} \
        --validateMappings -o $outdir >> {log} 2>&1
fi
if [ "$(basename {output.quant})" != "quant.sf" ]; then
    mv -f $outdir/quant.sf {output.quant}
fi
"""
