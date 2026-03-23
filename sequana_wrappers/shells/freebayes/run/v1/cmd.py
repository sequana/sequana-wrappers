CMD = """\
freebayes {params.options} --ploidy {params.ploidy} \
    -f {input.ref} -v {output.vcf} {input.bam} > {log} 2>&1
"""
