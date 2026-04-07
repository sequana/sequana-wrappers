CMD = """\
snpeff_dir=$(mktemp -d)
mkdir -p $snpeff_dir/data/{params.genome_name}
cp {input.ann} $snpeff_dir/data/{params.genome_name}/genes.gbk
echo "data.dir = $snpeff_dir/data" > $snpeff_dir/snpEff.config
echo "{params.genome_name}.genome : {params.genome_name}" >> $snpeff_dir/snpEff.config
snpEff build -c $snpeff_dir/snpEff.config -genbank -v {params.genome_name} > {log} 2>&1
snpEff ann {params.options} -c $snpeff_dir/snpEff.config \
    -csvStats {output.csv} \
    -stats {output.html} \
    {params.genome_name} {input.vcf} > {output.vcf} 2>> {log}
rm -rf $snpeff_dir
"""
