# Documentation

Sequana coverage detects and characterises automatically low and high 
genome coverage regions.
It provides a useful HTML report with dynamic plot and table.
Moreover, a CSV file with all metrics computed is created.
This CSV can be used to regenerate the sequana_coverage report.


**Required input:**

- **bed**: a BED file (built e.g with samtools_depth -aa input.bam)
- **fasta**: a FASTA file of the reference.
- **gbk**: a GENBANK file

**Required output:**

- FASTA file with locus names.

**Log:**

- a log file 

# Configuration

    ##############################################################################
    #
    # :Parameters:
    #
    # :param circular: is your genome circular or not ?
    # :param chunksize: for large genomes, split the data into chunks
    # :param double_threshold: double threshold for clustering. Keep 0.5 if you do
    #     not know. Otherwise, checkout the online documentation on
    #     sequana.readthedocs.io
    # :param high_threshold: keep 4 or check the online documentation
    # :param low_threshold: keep -4 or check the online documentation
    # :param mixture_models: keep to 2.
    # :param window: the W parameter of the running median. Keep as long as twice
    #     the deleted/depleted/duplicated you want to identify or to avoid. short
    #     genome will be set to genome length divided by 5 automatically. 
    sequana_coverage:
        do: true
        circular: true
        chunksize: 6000000
        double_threshold: 0.5
        gc_window_size: 201
        high_threshold: 4.0
        low_threshold: -4.0
        mixture_models: 2
        window_size: 3001

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

