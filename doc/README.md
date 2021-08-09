# The Sequana Wrapper Repository

The Sequana Wrapper Repository is a collection of reusable wrappers used in Sequana Pipelines inspired from [Snakemake Wrapper Repository](https://github.com/snakemake/snakemake-wrappers/).



# Notes for developers


- Rules may use the threads key, which can be defined in the general config file.
- the **params** section should contain an **options** also define in the config
  file.
- config file should prefix keys related to directories and files using the *_directory* or *_file* suffices


Consider this example:: 

    rule falco:
        input: whatever
        output: whatever
        log:
            "samples/{sample}/falco.log"
        threads:
            config['falco']['threads']
        params:
            options=config['falco']['options'],
            wkdir=config['falco']['working_directory']
        wrapper:
            "falco/wrappers/falco"


