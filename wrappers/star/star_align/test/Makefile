all:
	snakemake -s Snakefile --forceall -c 1 --wrapper-prefix git+file://$(realpath $(CURDIR)/../../../ )

clean:
	rm -rf .snakemake
	rm -rf star_index
	rm -rf star_mapping
	rm -rf star_mapping2
