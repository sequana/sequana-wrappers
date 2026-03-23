CMD = """\
outdir=$(dirname {output})
prefix=$(basename {output} .gff)
prokka --force {params.options} \
    --cpus {threads} \
    --outdir $outdir \
    --prefix $prefix \
    {input.assembly} > {log} 2>&1
"""
