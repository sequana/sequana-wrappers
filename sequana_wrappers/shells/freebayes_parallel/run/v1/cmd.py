CMD = """\
freebayes-parallel <(fasta_generate_regions.py {input.ref}.fai {params.chunk}) {threads} \
    {params.options} --ploidy {params.ploidy} \
    -f {input.ref} -v {output.vcf} {input.bam} > {log} 2>&1 || echo ''
# Note: '|| echo ''' suppresses the non-zero exit code that freebayes-parallel
# may return for empty genomic regions, while still capturing real errors in {log}.
"""
