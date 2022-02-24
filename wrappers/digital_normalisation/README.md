# Documentation

Digital normalisation is a method to normalise coverage of a sample in
fixed, low memory and without any reference.
The assembly with normalised data provides results qs good or even better than
assembling the unnormalised data.
Furthermore, SPAdes with normalised data is notably speeder and cost less
memory than without digital normalisation.


**Required input:**

- Fastq files (gzip compressed or not). Please provide a list for paired data, a string otherwise

**Required output:**

- Fastq files uncompressed. Please provide a list for paired data, a string otherwise

Parameters:

**Log:**

-  Log with stdout and sterr of Khmer


# Configuration


    ##############################################################################
    # Khmer - Digital Normalisation
    #
    # :Parameters:
    #
    # - do: if unchecked, this rule is ignored.
    # - ksize: kmer size used to normalised the coverage.
    # - cutoff: when the median k-mer coverage level is above this number the read
    #       is not kept.
    # - max_memory_usage: maximum amount of memory to use for data structure.
    # - threads: number of threads to be used.
    # - options: any options recognised by normalize-by-median.py.
    #
    digital_normalisation:
        do: yes
        ksize: 20
        cutoff: 20
        max_memory_usage: 4e9
        threads: 4
        options: ''

# Example

Note that the input and output must be named

:

    rule digital_normalisation_PE:
        input:
            fastq=["data_R1_.fastq.gz", "data_R2_.fastq.gz"]
        output:
            fastq=["data_R1_.dnpe.fastq", "data_R2_.dnpe.fastq"],
        log:
            "dn_PE.log"
        params:
            prefix = "dn",
            ksize = 20,
            cutoff = 20,
            m = 1000000000,
            options = '',
        threads: 1
        wrapper:
            "main/wrappers/digital_normalisation"


# References

