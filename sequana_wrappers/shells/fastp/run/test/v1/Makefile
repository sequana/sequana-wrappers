PARENT := $(shell realpath ..)

all:
	snakemake -s Snakefile -c4 --use-singularity --forceall \
		--singularity-args "-B $(PARENT)"

clean:
	rm -rf trimmed report logs .snakemake
