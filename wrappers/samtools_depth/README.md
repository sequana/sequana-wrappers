# Documentation

Samtools Depth creates a BED file with the coverage depth for each base
position. It can also compute multiple BAM files and concatenate results in
one BED file.

**Required input:**

- A Sorted BAM file or list of bam file.

**Required output:**
  
- A BED file with coverage for each base.

**Required parameters:**

- **options**: a list of valid Samtools Depth options

**Log:**

- The redirected standard error of Samtools Depth

# Configuration

	######################################################################
	# Samtools Depth section
	#
	# :Parameters:
	#
	# - options: string with any valid Samtools Depth options
    samtools_depth:
        options: ''

# Example

    rule samtools_depth:
        input:
            "{sample}/bamfile/{sample}.sorted.bam"
        output:
            "{sample}/samtools_depth/{sample}.bed" 
        log:
            "{sample}/samtools_depth/samtools_depth.log"
        params:
            options = config['samtools_depth']['options']
        wrapper:
            "main/wrappers/samtools_depth"

# References

- https://www.htslib.org/doc/samtools-depth.html
