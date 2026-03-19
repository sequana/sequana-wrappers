all:
	snakemake -s Snakefile -c4 --use-singularity --forceall

clean:
	rm -rf genome.fasta.bwt genome.fasta.fai genome.fasta.pac genome.fasta.ann \
	       genome.fasta.amb genome.fasta.sa logs .snakemake
