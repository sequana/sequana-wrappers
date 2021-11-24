# Documentation

Sequana coverage detects and characterises automatically low and high 
genome coverage regions.
It provides a useful HTML report with dynamic plot and table.
Moreover, a CSV file with all metrics computed is created.
This CSV can be used to regenerate the sequana_coverage report.


Required input:

- BED file (built e.g with samtools_depth -aa input.bam)
- FASTA file of the reference.
- GENBANK file

Required output:

- FASTA file with locus names.

Log:

- a log file 

# Configuration

    sequana_coverage:
        do: true
        circular: true
        window_size: 3001
        chunksize: 5000000
        double_threshold: 0.5
        gc_window_size: 201
        high_threshold: 4.0
        low_threshold: -4.0
        mixture_models: 2

# Example

    rule sequana_coverage:
        input:
            bed = "test.bed",
            fasta = "measles.fa",
            gbk = "measles.gbk"
        output:
            html = "test/sequana_coverage.html"
        params:
            mixture_models = 2,
            window_size = 3001,
            double_threshold = 0.5,
            circular = True,
            chunksize = 5000000,
            gc_size = 201,
        wrapper:
            "main/wrappers/sequana_coverage"

