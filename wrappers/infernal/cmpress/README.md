# Documentation

This wrapper calls **cmpress** to compress a covariance model (CM) database
for use with **cmscan**. It constructs binary compressed datafiles from a CM
file in standard Infernal-1.1 format.

**Required input:**

- A covariance model file (`.cm`)

**Required output:**

- Four compressed binary index files: `.cm.i1i`, `.cm.i1f`, `.cm.i1m`, `.cm.i1p`

**Optional parameters:**

- **extra**: any additional valid cmpress options (default to "")

**Log:**

- a log file capturing stderr from cmpress

# Configuration

    ##############################################################################
    # Infernal cmpress section.
    #
    # :Parameters:
    #
    # - extra: string with any additional valid cmpress options
    #
    infernal_cmpress:
        extra: ""

# Example

    rule infernal_cmpress:
        input:
            "database.cm"
        output:
            "database.cm.i1i",
            "database.cm.i1f",
            "database.cm.i1m",
            "database.cm.i1p"
        log:
            "logs/cmpress.log"
        params:
            extra=""
        wrapper:
            "main/wrappers/infernal/cmpress"

# References

- http://eddylab.org/infernal/
- Nawrocki EP and Eddy SR. Infernal 1.1: 100-fold faster RNA homology searches.
  Bioinformatics 29:2933-2935 (2013).
