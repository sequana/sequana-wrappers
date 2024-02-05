# Documentation

Same as the freebayes wrapper except that it is a parallelised version.
Only 1 additional option called **chunk** that is set to 1,000,000 by default. 
And the **threads** key.

# Configuration

	######################################################################
	# Freebayes section
	#
	# :Parameters:
	#
	# - ploidy: sets the default ploidy for the analysis
	# - options: string with any valid Freebayes options
    freebayes:
		ploidy: 1
		options: "--legacy-gls"
        threads: 4

# Example

    rule freebayes:
        input:
            bam = "{sample}/bamfile/{sample}.sorted.bam",
            ref = "measles.fa"
        output:
            vcf = "{sample}/freebayes/{sample}.vcf"
        log:
            "{sample}/freebayes/freebayes.log"
        params:
            ploidy = config["freebayes"]["ploidy"],
            options = config["freebayes"]["options"]
            chunk = 5000
        threads: 
            4
        wrappers"
            "main/wrappers/freebayes_parallel/

# References

- https://github.com/freebayes/freebayes
