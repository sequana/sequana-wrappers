# Documentation

Freebayes is a variant caller designed to find SNPs and short INDELs from a
BAM file. It produces a very well-annotated VCF output.
Moreover, it provides a quality score calculated by a bayesian model.


Required input:

- Sorted BAM file.
- FASTA file of the reference genome.

Required output:

- VCF file of detected variants.

Log:

- a log file

# Configuration


    freebayes:
       ploidy: 1
       options: --legacy-gls


# Example

    rule freebayes:
        input:
            bam = "test.sorted.bam"
            ref = "measles.fa"
        output:
            vcf = "test.vcf"
        log:
            "freebayes.log"
        params:
            ploidy = 1,
            options = ""
        wrappers"
            "main/wrappers/freebayes/
