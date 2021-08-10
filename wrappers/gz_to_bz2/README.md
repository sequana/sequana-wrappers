rule gz_to_bz2:
    """**Convert fastq.gz files to fastq.bz2 files**

    Here are the steps followed by the rule. Any failure stops the
    process and the original file is untouched. If all succeed, the
    input is deleted.

        #. the input GZ file is checked for integrity.
        #. the input GZ file is decompressed with **pigz** and redirected
           a pipe to **pbzip2** executable into a BZ2 output.
        #. the output is checked for integrity with **pbzip2**.
        #. the input GZ file is deleted.

    Required input:
         - the input gzipped files (wildcards possible)
    Required output:
         - the output bzipped files (wildcards possible)
    Required parameters:
        - threads:


    configuration requirements::

        compressor:
            - threads

    """
    input: "{dataset}.gz"
    output: "{dataset}.bz2"
    threads: config['compressor']['threads']
    wrapper: main/wrappers/gz_to_bz2
