# Documentation

Sequana taxonomy performs taxonomic analysis using Kraken2. 
It is essentially a wrapper of kraken that hides technical details
and a single standalone to cope with paired/unpaird data, multiple
kraken databases, conversion of taxons, etc. It also provides a HTML
report

**Required input:**

- **fastq**: a BED file (built e.g with samtools_depth -aa input.bam)

**Required output:**

The output of sequana_taxonomy is multiple. In the output directory, 
you will find a summary.html, and a sub directory called kraken with
a CSV file called kraken.csv, a summary file called kraken.out.summary
a json summary file called summary.jsona and finally unclassified 
fastq files. The option store_unclassified is currently required. If
set to False, empty files are created.



# Configuration

    ##############################################################################
    #
    # :Parameters:
    #
    #
    # * paired indicates whether input fastq is paired or not.
    # * databases: a list of valid databases
    # * store_unclassified
    #
    # optional ones:  confidence (default to 0), level (default to INFO) and
    #                 options (default to "")
    sequana_taxonomy:
        databases:
            - /home/user/.config/sequana/kraken2_dbs/viruses_masking
        level: INFO
        confidence: 0.05
        store_unclassified: false
        options: ''
        threads: 4

# Example

    rule sequana_taxonomy:
        input:
            "test_R1_.fastq"  # second file for paired data is possible
        output:
            html         = '{sample}/summary.html',
            csv          = '{sample}/kraken/kraken.csv',
            summary      = '{sample}/kraken/kraken.out.summary',
            summary_json = '{sample}/kraken/summary.json',
            unclassified = '{sample}/kraken/unclassified.fastq'
        threads:
            4
        params:
            # required
            paired=True,
            databases=['toydb'],  # a list of valid kraken databases
            store_unclassified=True,
            # optional
            confidence=0,
            level="INFO",
            options=""
        container:
            https://zenodo.org/record/7963917/files/sequana_tools_0.15.1.img
        wrapper:
            "main/wrappers/sequana_taxonomy"





