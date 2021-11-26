# Documentation

dsrc_to_gz converts fastq.dsrc files to fastq.gz files

Here are the steps followed by the rule. Any failure stops the
process and the original file is untouched. If all succeed, the
input is deleted.

- the input DSRC file is decompressed with **dsrc** and redirected a pipe to
  **pigz** executable into a GZ output.
- the output is checked for integrity with **pigz**.
- the input DSRC file is deleted.

**Required input:**

- a dsrc compressed FASTQ file

**Required output:**

- a gz compressed FASTQ file

# Configuration

	######################################################################
	# dsrc_to_gz section
	#
	# :Parameters:
	#
	dsrc_to_gz:
		threads: 4
		options: "-m2"

# Example

	rule dsrc_to_gz:
		input:
			"{sample}.fq.dsrc"
		output:
			"{sample}/dsrc_to_gz/{sample}.fq.gz"
		params:
			options = config["dsrc_to_gz"]["options"]
		threads:
			config["dsrc_to_gz"]["threads"]
		wrapper:
			"main/wrappers/dsrc_to_gz"
