# Documentation


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

    rule add_read_group:
        input:
            fasta="test.fasta",
            ann="test.gbk"
        output:
            "{sample}/bamtools/{sample}.sorted.bai
        log:
            "common_logs/snpeff_add_locus_in_fasta.log"
        wrapper:
            "main/wrappers/snpeff_add_locus_in_fasta
