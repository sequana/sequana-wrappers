# Documentation

This wrapper calls **cmscan** to search sequences against collections of
covariance models (CMs) that have been prepared with **cmpress**. It is used
to identify RNA structural homologs in sequence data.

**Required input:**

- **fasta**: the input sequence file to search
- **profile**: one of the binary index files produced by cmpress (e.g. `.cm.i1i`)

**Optional output:**

- **outfile**: the main human-readable output file. If not specified, output goes to `/dev/null`.
- **tblout**: a tab-separated table of per-target hits.

**Required parameters:**

- **evalue_threshold**: report hits with E-value ≤ this value (default: 10).
  Ignored if `score_threshold` is set.
- **score_threshold**: report hits with bit score ≥ this value. Overrides
  `evalue_threshold` when set.
- **extra**: any additional valid cmscan options (default to "")

**Log:**

- a log file capturing stderr from cmscan

# Configuration

    ##############################################################################
    # Infernal cmscan section.
    #
    # :Parameters:
    #
    # - evalue_threshold: E-value threshold for reporting hits
    # - extra: string with any additional valid cmscan options
    #
    infernal_cmscan:
        evalue_threshold: 10
        extra: ""
        threads: 4

# Example

    rule infernal_cmscan:
        input:
            fasta="transcripts.fa",
            profile="rfam.cm.i1i"
        output:
            tblout="rfam-hits-tblout.txt"
        log:
            "logs/cmscan.log"
        params:
            evalue_threshold=10,
            extra=""
        threads: 4
        wrapper:
            "main/wrappers/infernal/cmscan"

# References

- http://eddylab.org/infernal/
- Nawrocki EP and Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches.
  Bioinformatics 29:2933-2935 (2013).
