# Documentation

gz_to_bz2 Converts fastq.gz files to fastq.bz2 files.

Here are the steps followed by the rule. Any failure stops the
process and the original file is untouched. If all succeed, the
input is deleted.

- the input GZ file is checked for integrity.
- the input GZ file is decompressed with **pigz** and redirected a pipe to
  **pbzip2** executable into a BZ2 output.
- the output is checked for integrity with **pbzip2**.
- the input GZ file is deleted.

**Required input:**

- A FASTQ gzip compressed file

**Required output:**

- A FASTQ bz2 compressed file

# Configuration

	######################################################################
	# gz_to_bz2
	#
	gz_to_bz2:
		threads: 4

# Example

	rule gz_to_bz2:
		input:
			"{sample}.fq.gz"
		output:
			"{sample}/gz_to_bz2/{sample}.fq.bz"
		threads:
			config["gz_to_bz2"]["threads"]
		wrapper:
			"main/wrappers/gz_to_bz2"
