rule fastqc:
    """Calls FastQC on input data sets (paired or not)

    This wrapper takes 1 or 2 FastQ files as input. It creates a **fastqc.done**
    file defined by the user. Result are stored in the directory extracted from.
    the output file name. In other words, if the output is defined as
    fastqc/fastqc.done, results will be stored in ./fastqc directory. 

    Required input:
        - one file or two files (if paired)

    Required output:
        - a filename used as a trigger of job done. Be aware that the root
          directory of this file is used to store the results

    Required parameters
        - options: a list of valid fastqc options
        - working_directory: place where results of fastqc are saved

    Log:
        - a log file is created

    Required configuration section:
        .. code-block:: yaml

            fastqc:
                options: "-nogroup"   # a string with fastqc options
                threads:              # optional. if not set, 4 threads are used

    References:
        - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """
    input:
        "{sample}_R1_.fastq", "{sample}_R2_.fastq"
    output:
        done = "samples"/{sample}/fastqc.done"
    params:
        options = config['fastqc'][options']
        working_directory = "samples/{sample}"
    threads: 
        config['fastqc']["threads"]
    log:
        "fastqc/fastqc.log"
    wrapper:
        "main/wrappers/fastqc"


