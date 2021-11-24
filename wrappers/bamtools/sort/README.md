# Documentation

Sorting of a BAM file using bamtools

**Required input:**

- the input BAM file

**Required output:**

- the sorted BAM file

**Log:**

- a log file produced by bamtools

# Configuration

There is no configuration required for this wrapper.
However, to be consistent with other wrapper, the
params.options is used (empty by default).

# Example

    rule bamtools_index:
        input:
            "{sample}/bamtools/{sample}.bam
        output:
            "{sample}/bamtools/{sample}.sorted.bam
        log:
            "{sample}/bamtools/bamtools_sort.log"
        wrapper:
            "main/wrappers/bamtools/sort"

