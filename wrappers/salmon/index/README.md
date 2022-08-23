# Documentation

This wrapper calls salmon to build index of a genome.
Output is saved in a working directory infered from the 
parent directory of the output file. If this is not what you want,
you can set the **wkdir** parameter

**Required input:**

- **fasta**: genome in FastA format
- **gff**: annotation in GFF format 

**Required output:**

- **done**: a filename used as a trigger the end of the processing.

Since there is only one output, you may name it differently

**Optional parameters:**

- **options**: a list of valid salmon index options (Default to "")
- **threads**: thread number (default to 1)
- **wkdir**: default to the parent directory o the output.done file

**Log:**

- a log file common to gffread and salmon index

# Configuration

    ##############################################################################
    # Salmon
    #
    # :Parameters:
    #
    # - options: string with any valid Salmon index options
    salmon_index:
        options: ""
        threads: 2

# Example

without annotation::

    rule salmon_index:
        input:
            fasta = "genome.fasta"
            gff = "genome.gff"
        output:
            done = "salmon/salmon.done"
        params:
            options="",
        log:
            "salmon/salmon.log"
        wrapper:
            "main/wrappers/salmon_index"


# References

