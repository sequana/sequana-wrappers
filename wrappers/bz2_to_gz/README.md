rule bz2_to_gz:
    """**Convert fastq.gz files to fastq.bz2 files**

    Here are the steps followed by the rule. Any failure stops the
    process and the original file is untouched. If all succeed, the
    input is deleted.

        #. the input BZ2 file is checked for integrity.
        #. the input BZ2 file is decompressed with **pbunzip2** and redirected
           a pipe to **pigz** executable into a GZ output.
        #. the output is checked for integrity with **pigz**.
        #. the input BZ2 file is deleted.

    Required input:
         - the input bzipped files (wildcards possible)
    Required output:
         - the output gzipped files (wildcards possible)
    Required parameters:
        - threads:


    configuration requirements::

        compressor:
            - threads

    """
    input: "{dataset}.bz2"
    output: "{dataset}.gz"
    threads: config['compressor']['threads']
    wrapper: main/wrappers/bz2_to_gz
