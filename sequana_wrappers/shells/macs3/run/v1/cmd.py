CMD = """\
format_opt="-f BAM"
if [ "{params.paired}" = "True" ]; then format_opt="-f BAMPE"; fi
broad_opt=""
if [ "{params.mode}" = "broad" ]; then
    broad_opt="--broad --broad-cutoff {params.broad_cutoff}"
fi
macs3 callpeak -B --SPMR \
    -t {input.inputs} \
    -c {input.controls} \
    -g {params.genome_size} \
    -n {params.prefix} \
    --bw {params.bandwidth} \
    {params.options} \
    -q {params.qvalue} \
    $format_opt \
    $broad_opt \
    --outdir {params.outdir} > {log} 2>&1
"""
