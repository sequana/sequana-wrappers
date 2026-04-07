CMD = """\
snpeff_dir=$(mktemp -d)
config_file=$(mktemp)
echo "# SnpEff config file" > $config_file
echo "{params.genome_name}.genome : {params.genome_name}" >> $config_file
mkdir -p $snpeff_dir/{params.genome_name}
cp {input.ann} $snpeff_dir/{params.genome_name}/genes.gbk
snpEff build -c $config_file -dataDir $snpeff_dir -genbank -noCheckProtein -v {params.genome_name} > {log} 2>&1
snpEff ann {params.options} -c $config_file -dataDir $snpeff_dir \
    -csvStats {output.csv} \
    -stats {output.html} \
    {params.genome_name} {input.vcf} > {output.vcf} 2>> {log}
rm -f $config_file
rm -rf $snpeff_dir
"""
