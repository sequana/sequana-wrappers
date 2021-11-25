# Documentation

This wrapper makes long reads correction and/or assembly using Canu.

This wrapper takes fastq or fasta file. Using `corrected` or `trimmed` step option,
the output is corrected reads in a fasta.gz file. Without step option, the output
is the assembly in a fasta file.

**Required input:**

- Long reads file (FASTQ/FASTA format)

**Required output:**

- Assembly result in fasta or corrected reads in fasta.gz.
- The `canu{step}.done` file is a trigger file to let Snakemake knows that canu computation is over.
 
**Required parameters**

- **preset**: Any preset in this list:
  - `"pacbio"`: Raw pacbio data 
  - `"pacbio-hifi"`: Hifi pacbio data
  - `"nanopore"`: Nanopore data
- **genome_size**: The expected genome size.
- **step**: The step that you want to do in this list:
  - `"-correct"`: Canu read correction.
  - `"-trim"`: Canu read trimming.
  - `""`: Default Canu assembly.
- **use_grid**: Let Canu handle the usage of your cluster or not.
- **options**: a list of valid Canu options.


# Configuration

    ##############################################################################
    # Canu long read assembly
    #
    # :Parameters:
    #   
    # - preset: Any preset in this list. (pacbio, pacbio-hifi, nanopore)
    # - genome_size: An estimate of the size of the genome. Common suffices are allowed, for example, 3.7m or 2.8g.
    # - use_grid: let canu run steps on cluster.
    # - step: Any step in this list. (-correct, -trim, "")
    # - options: any options recognised by canu
    # - threads: Number of threads to use
    canu:
        preset: 'pacbio'
        genome_size: '3.3m'
        step: ''
        use_grid: true
        options: ''
        threads: 1

# Example

    rule canu:
        input:
            "nice_pb_long_reads.fastq"
        output:
            "canu/nice_assembly.contigs.fasta",
            "canu/canu.done"
        params:
            preset = "pacbio",
            genome_size = "3G",
            step = "",
            use_grid = True,
            options = ""
        threads: 1
        wrapper:
            "main/wrappers/canu"

# References
- https://canu.readthedocs.io/en/latest/quick-start.html
