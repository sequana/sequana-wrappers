# Documentation

SnpEff adds annotation of variants detected in a VCF file. It annotates
using the old 'EFF' field instead of 'ANN' field. The latter does not
provide the codon change information.

**Required input:**

- **vcf**: VCF file of detected variants.
- **ann**: Annotation genbank file

**Required output:**

- **vcf**: Annotated VCF file.
- **html**: HTML report
- **csv**: CSV file with Variants

# Configuration

	######################################################################
	# SNPEff section
	#
	# :Parameters:
	#
	# - options: string with any valid SNPEff options
	snpeff:
		options: '-no-downstream'

# Reference

- https://pcingola.github.io/SnpEff/se_introduction/


# Example

    rule snpeff:
        input:
            vcf = "{sample}/freebayes/{sample}.vcf",
            ann = "measles.gbk"
        output:
            vcf = "{sample}/snpeff/{sample}.vcf",
            csv = "{sample}/snpeff/{sample}.csv",
            html = "{sample}/snpeff/{sample}.html"
        log:
            "{sample}/snpeff/{sample}.log"
        params:
            options = config["snpeff"]["options"]
        wrapper:
            "main/wrappers/snpeff"

