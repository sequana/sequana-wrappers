# Documentation

SnpEff requires the locus names in the annotation file and in the FASTA
file (contig name) to be identical. To make this is true, this rule adds
locus names of the genbank file into the FASTA file before the mapping.

**Required input:**

- **fasta** FASTA file of the reference.
- **ann** GENBANK or GFF file for annotation

**Required output:**

- FASTA file with locus names.

**Log:**

- a log file 

# Configuration

There is no configuration required for this wrapper.

# Example

    rule snpeff_add_locus_in_fasta:
        input:
            fasta="{sample}.fas",
            ann="{sample}.gbk"
        output:
            "{sample}/snpeff_add_locus_in_fasta/{sample}.fas
        log:
            "{sample}/snpeff_add_locus_in_fasta/{sample}.log"
        container:
            "https://zenodo.org/record/7963917/files/sequana_tools_0.15.1.img"
        wrapper:
            "main/wrappers/snpeff_add_locus_in_fasta
