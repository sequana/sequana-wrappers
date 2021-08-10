rule falco:
    """Calls Faslco on input data sets (paired or not)

    This wrapper takes 1 or 2 FastQ files as input and runs Falco on the files. 

    Required input:
        - one file or two files (if paired)

    Required output:
        - one output file from falco such as summary.txt 

    Required parameters
        - options: a list of valid falco options (see below)
        - working_directory: place where results will be stored
    
    Threads:
        - threads are taken from the config file (see below)

    Log:
        - a log file is created

    Required configuration section:
        .. code-block:: yaml

            falco:
                options: " "   # a string with falco options
                threads:       # can be used if set

    References:
        - http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    """
    input:
        "{sample}_R1_.fastq", "{sample}_R2_.fastq"
    output:
        done = "samples"/{sample}/summary.txt"
    params:
        options = config['falco'][options']
        working_directory = "samples/{sample}"
    threads: 
        config['falco']["threads"]
    log:
        "logs/falco.log"
    wrapper:
        "main/wrappers/falco"


