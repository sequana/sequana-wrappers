# Documentation

This wrapper creates images and statistics related to FastQ data

This wrapper takes 1 FastQ file as input.

**Required input:**

- one FastQ file.

**Required output:**

- **gc**:  PNG image of the GC content
- **json**:  some stats such as GC content, mean read length, etc
- **boxplot**:  a fastqc-like quality image


**Required parameters:**

- **max_reads**: uses only 500,000 reads

# Configuration

Not required

# Example

    rule fastq_stats:
        input:
            "{sample}.fastq.gz"
        output:
            gc="{sample}/fastq_stats/{sample}_gc.png",
            boxplot="{sample}/fastq_stats/{sample}_boxplot.png",
            json="{sample}/fastq_stats/{sample}.json"
        params:
            max_reads=500000
        wrapper:
           "main/wrappers/fastq_stats"

# References

- https://sequana.readthedocs.io
