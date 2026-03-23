CMD = """\
snpeff_dir=$(mktemp -d)
snpEff build -dataDir $snpeff_dir -genbank -v {params.genome_name} {input.ann} > {log} 2>&1
snpEff ann {params.options} -dataDir $snpeff_dir \
    -csvStats {output.csv} \
    -stats {output.html} \
    {params.genome_name} {input.vcf} > {output.vcf} 2>> {log}
rm -rf $snpeff_dir
"""
