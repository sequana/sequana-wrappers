# The Sequana Wrapper Repository

The Sequana Wrapper Repository is a collection of reusable wrappers used in the [Sequana Pipelines](https://sequana.readthedocs.io). This repository is inspired from [Snakemake Wrapper Repository](https://github.com/snakemake/snakemake-wrappers/). All the Sequana wrappers can can be used in your Snakemake pipelines as well.


[![Tests](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml/badge.svg)](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml)
[![Tests](http://joss.theoj.org/papers/10.21105/joss.00352/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00352)


|||
| --- | --- |
| Overview | A set of Snakemake wappers to help building Sequana pipelines |
|Status | Production |
|Issues | Please fill a report on [github/sequana/sequana-wrappers](https://github.com/sequana/sequana/issues) |
|Python version | Python 3.7, 3.8, 3.9|
|Citation| Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16) ), 352,  [JOSS DOI doi:10.21105/joss.00352 ](http://www.doi2bib.org/bib/10.21105%2Fjoss.00352) |


# Usage

    snakemake --wrapper-prefix https://github.com/sequana/sequana-wrappers

or with a local copy in /home/user:

    cd
    git clone git@github.com:sequana/sequana-wrappers.git sequana_wrappers
    snakemake --wrapper-prefix git+file:////home/user/sequana_wrappers


# Notes for developers

## Overview

The wrappers directory contains the wrappers. Each sub-directory is dedicated to
a wrapper related to a given software/application. 

Here is an example of a wrapper tree structure:

    fastqc
    ├── environment.yaml
    ├── README.md
    ├── test
    │   ├── README.md
    │   ├── Snakefile
    │   ├── test_R1_.fastq
    │   └── test_R2_.fastq
    └── wrapper.py

Note that some software may have several sub wrappers (see the bowtie1 wrapper for instance).

A wrapper directory must contain a file called **wrapper.py** where the
developers must provide the core of the wrapper. There is no specific
instructions here except to write good code as much as possible (with comments).

A wrapper directory should have a **test** directory for continuous integration
with a **Snakefile** to be tested and possibly data file **Do not add large files here**. A
**README.md** should be added to explain the origin of the test data files.
Finally, include your tests in the main [**test.py**](test.py) file
of the root of the repository (not the wrapper itself).

For testing purposes, you should also add a file called **environment.yaml**
to tell what are the required packages to be installed for the test (and wrapper)
to work.

Finally, for the documentation, we ask the developer to create a **README.md** file
described here below.

To test your new wrapper (called *example* here), type:

   pytest test.py -k test_example

## The config file

If required in a wrapper, parameters must be defined in a **config.yaml** file.
Similarly for threading. Consider the following pointswhen writting a wrapper:

- The thread paramter should also be a parameter in config file.
- the **params** section should contain a key called **options** also define in the config file.
- keys or parameters related to directories and files should use the *_directory* or *_file* suffices. This is for
  Sequanix application to automatically recognised those options with a dedicated widget.

Consider this example:

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

You config file will look like:

    falco:
        threads: 4
        options="--verbose"
        working_directory: 'here'


## Naming arguments of the different sections

In all sections (e.g., input, params), if there is only one input, no need to name it, otherwise, please do.

    rule example1:
        input:
            "test.bam"
        output:
            "test.sorted.bam"
        ...

but:

    rule example1:
        input:
            "test.bam"
        output:
            bam="test.sorted.bam"
            bai="test.sorted.bam.bai"
        ...




## Documentation

Each wrapper should have a dedicated documentation explaining the input/output with a usage example. It should also document the expected configuration file.  The file must be formatted in markdown. It must contain a **Documentation** and **Example** sub sections. If a **Configuration** section is found, it is also added to the documentation. This **README.md**  file will be rendered automatically via a Sequana sphinx plugin. Consider the [fastqc](wrappers/fastqc/README.md) directory for a workable example rendered [here](https://sequana.readthedocs.io/en/master/wrappers.html#fastqc).



## Faqs

When adding new recipe / testing, you may face several issues.

Make sure you have added/commited the files you want to test.

Generally, you should create a branch, add the recipe. In the test file, when setting the wrapper
you should replace the first 'main' by your branch name and then you can test the wrapper locally as
follows:

   cd wrapper/your_recipe/test
   snakemake -s Snakefile  -j 1 --wrapper-prefix git+file:///YOURPATH/sequana-wrappers/ -f -p

Once it works, do not forget to replace the branch name in the test.Snakefile with "main"





