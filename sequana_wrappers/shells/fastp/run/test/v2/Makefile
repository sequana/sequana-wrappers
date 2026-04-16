PARENT := $(shell realpath ..)

all:
	snakemake -s Snakefile -c4 --use-singularity --forceall \
		--singularity-args "-B $(PARENT)"

clean:
	rm -rf trimmed_pe trimmed_se report_pe report_se logs .snakemake
