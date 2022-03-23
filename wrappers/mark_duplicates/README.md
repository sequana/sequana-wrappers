# Documentation

This wrapper marks duplicates using picard tools so that tools such as variant
caller are aware of duplicated reads

**Required input:**

- BAM file

**Required output:**

- The first output must be the expected output BAM file with the 
  marked duplicated reads.
- The second output contains metrics

**Log:**

- a log file with stdout/stderr

**Required parameters:**

- remove_dup if set to true will remove the duplicated reads (default is False)
- tmpdir is used to stored temporary files


# Configuration


    #############################################################################
    # mark_duplicates (picard-tools) allows to mark PCR duplicate in BAM files
    #
    # :Parameters:
    #
    # - do: if unchecked, this rule is ignored. Mandatory for RNA-SeQC tool.
    # - remove: If true do not write duplicates to the output file instead of writing them with
    #            appropriate flags set.  Default value: false. This option can be set to 'null' to clear
    #            the default value. Possible values: {true, false}
    # - tmpdir: write tempory file on this directory (default TMP_DIR=/tmp/, but could be "TMP_DIR=/local/scratch/")
    #
    mark_duplicates:
        do: false
        remove: false ## may be True
        tmpdir: ./tmp/

# Example


    rule mark_duplicates:
        input:
            "test.bam"
        output:
            bam = "test/md.bam",
            metrics = "test/md.metrics"
        log:
            out = "test/log.out",
            err = "test/log.err"
        params:
            remove_dup = "false",
            tmpdir = "test/tmp"
        wrapper:
            "main/wrappers/mark_duplicates"


