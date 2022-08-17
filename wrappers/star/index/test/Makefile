all:
	snakemake -s Snakefile --forceall -c 1 --wrapper-prefix git+file://$(realpath $(CURDIR)/../../../ )

clean:
	rm -rf .snakemake
	rm -rf star_index
	rm -rf star_index2
	rm -rf star_index3
