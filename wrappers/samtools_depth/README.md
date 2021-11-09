# Documentation

Samtools Depth creates a BED file with the coverage depth for each base
position. It can also compute multiple BAM files and concatenate results in
one BED file.


Required input:

- A Sorted BAM file or list of bam file.

- Required output:
  
- A BED file with coverage for each base.

Log:

- a log file with stdout/stderr

If provided, options field in the *params* sections is used.


# Configuration

    samtools_depth:
        options:


# Example


    rule samtools_depth:
        input:
            "test.sorted.bam"
        output:
            "test.bed" 
        log:
            "out.log"
        params:
            options=config['samtools_depth']['options']
        wrapper:
            "main/wrappers/samtools_depth"




