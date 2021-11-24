# Documentation

bz2_to_gz converts fastq.gz files to fastq.bz2 files

Here are the steps followed by the rule. Any failure stops the
process and the original file is untouched. If all succeed, the
input is deleted.

- the input BZ2 file is checked for integrity.
- the input BZ2 file is decompressed with **pbunzip2** and redirected a pipe to
  **pigz** executable into a GZ output.
- the output is checked for integrity with **pigz**.
- the input BZ2 file is deleted.

**Required input:**

- bzipped files (wildcards possible)

**Required output:**

- output gzipped files (wildcards possible)

# Configuration

	######################################################################
	# bz2_to_gz section
	#
	# :Parameters:
	#
	bz2_to_gz:
		threads: 4

# Example

	rule bz2_to_gz:
		input:
			"{sample}.fq.bz2"
		output:
			"{sample}/bz2_to_gz/{sample}.fq.gz"		
		threads:
			config["bz2_to_gz"]["threads"]
		wrapper:
			"main/wrappers/bz2_to_gz"

