# Documentation

SnpEff requires the locus names in the annotation file and in the FASTA
file (contig name) to be identical. To make this is true, this rule adds
locus names of the genbank file into the FASTA file before the mapping.

Required input:

- FASTA file of the reference.
- GENBANK file

Required output:

- FASTA file with locus names.

Log:

- a log file 

# Configuration

snpeff:
    annotation_file:  # the genbank file
    options:    # result filters options

# Example

    rule snpeff_add_locus_in_fasta:
        input:
            fasta="test.fasta",
            ann="test.gbk"
        output:
            "{sample}/bamtools/{sample}.sorted.bai
        log:
            "common_logs/snpeff_add_locus_in_fasta.log"
        wrapper:
            "main/wrappers/snpeff_add_locus_in_fasta
