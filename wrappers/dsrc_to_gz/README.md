rule dsrc_to_gz:
    """**Convert fastq.dsrc files to fastq.gz files**

    Here are the steps followed by the rule. Any failure stops the
    process and the original file is untouched. If all succeed, the
    input is deleted.

        #. the input DSRC file is decompressed with **dsrc** and redirected
           a pipe to **pigz** executable into a GZ output.
        #. Third, the output is checked for integrity with **pigz**.
        #. Fourth, the input DSRC file is deleted.

    :param input: {dataset}.dsrc
    :param output: {dataset}.gz 
    :third-party executables: dsrc and pigz

    configuration requirements::

        compressor:
            - threads

    """
    input: "{dataset}.dsrc"
    output: "{dataset}.gz"
    params: "-m2"
    threads: config["compressor"]["threads"]
    wrappers: main/wrapper/dsrc_to_gz

