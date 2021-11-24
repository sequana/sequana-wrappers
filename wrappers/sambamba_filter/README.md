# Documentation

This rule uses sambamba view to filter reads with a mapping quality lower
than a threshold. It also removes reads with multiple occurrence.

**Required input:**

- BAM file

**Required output:**

- the filtered output BAM file.

Parameters:

- **threshold**: Quality threshold used for filtering
- **options**: a list of valid sambamba filter options

**Log:**

- two log files named .out and .err

# Configuration

	######################################################################
	# Sambamba filter section
	#
	# :Parameters:
	#
	# - options: string with any valid Sambamba filter options
	# - threshold: quality threshold used for filtering
    sambamba_filter:
        threshold: 30 # Mapping quality score threshold
        options: ''

# Example

    rule sambamba_filter:
        input:
            "{sample}/bamfile/{sample}.sorted.bam"
        output:
            "{sample}/sambamba_filter/{sample}.sorted.bam"
        log:
            out = "{sample}/sambamba_filter/log.out",
            err = "{sample}/sambamba_filter/log.err"
        params:
            threshold = config["sambamba_filter"]["threshold"],
            options = config["sambamba_filter"]["options"]
        wrapper:
            "main/wrappers/sambamba_filter"

# References

- https://github.com/biod/sambamba
