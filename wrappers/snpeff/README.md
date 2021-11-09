# Documentation

SnpEff adds annotation of variants detected in a VCF file. It annotates
using the old 'EFF' field instead of 'ANN' field. The latter does not
provide the codon change information.

Required input:

- VCF file of detected variants.
- Annotation genbank file

Required output:

- Annotated VCF file.
- HTML report
- CSV file with Variants


# Configuration

    snpeff:
        reference: 'genes.gb' # The genbank file with annotation of the reference.
        options: '-no-downstream' # Any options

Results filter options:

- *-no-downstream*: Do not show DOWNSTREAM changes
- *-no-intergenic*: Do not show INTERGENIC changes
- *-no-intron*:     Do not show INTRON changes
- *-no-upstream*:  Do not show UPSTREAM changes
- *-no-utr*:       Do not show 5_PRIME_UTR or 3_PRIME_UTR changes



# Example

    rule snpeff:
        input:
            bam = "test.raw.vcf",
            ann = "measles.gbk"
        output:
            vcf = "test.ann.vcf",
            csv = "test.ann.csv",
            html = "snpeff.html"
        log:
            "snpeff.log"
        params:
            options = ""
        wrapper:
            "main/wrappers/snpeff"

# Note

.. seealso:: snpeff_add_locus_in_fasta


