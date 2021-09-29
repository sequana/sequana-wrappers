# Documentation

Indexing of a BAM file using bamtools

Required input:

- the input BAM file

Required output:

- the indexed BAM file (same as input + .bai extension)

Log:

- a log file

# Configuration

There is no configuration required for this wrapper.
However, to be consistent with other wrapper, the 
params.options is used (empty by default).


# Example

    rule bamtools_index:
        input:
            "{sample}/bamtools/{sample}.sorted.bam
        output:
            "{sample}/bamtools/{sample}.sorted.bai
        log:
            "{sample}/bamtools/bamtools_index.log"
        wrapper:
            "main/wrappers/bamtools/index"

