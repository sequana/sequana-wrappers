# Documentation

Get stats about a BAM file using bamtools

**Required input:**

- the input BAM file

**Required output:**

- the output text file


# Configuration

There is no configuration required for this wrapper.
However, to be consistent with other wrappers, the
params.options may be used (empty by default).

# Example

    rule bamtools_stats:
        input:
            "{sample}/bamtools/{sample}.bam
        output:
            "{sample}/bamtools/{sample}.txt
        container:
            "https://zenodo.org/record/7345682/files/bamtools_2.5.2.img"
        wrapper:
            "main/wrappers/bamtools/stats"

