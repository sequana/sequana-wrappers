# Documentation

This wrapper annotates assembly with **Prokka**.

**Required input:**

- **assembly**: the input assembly fasta file.

**Required output::**

- **output**: any file you need (.gbk/.gff/.fna).
 
**Required parameters:**

- **options**: a list of valid prokka options.

**Log:**

- a log file generated by prokka is created.


# Configuration


    ##############################################################################
    # Prokka
    #
    # :Parameters:
    #
    # - options: any options recognised by prokka cli.
    # - threads: number of threads to be used.
    #
    prokka:
        options:
        threads: 4


# Example

    rule prokka:
        input: "assembly/contigs.fasta"
        output: "prokka/contigs.gff"
        params:
            options=config["prokka"]["options"]
        log:
            "logs/prokka.log"
        threads:
            config["prokka"]["threads"]
        wrapper:
            "main/wrappers/prokka"