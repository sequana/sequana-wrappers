# The Sequana Wrapper Repository

The Sequana Wrapper Repository is a collection of reusable wrappers used in the [Sequana Pipelines](https://sequana.readthedocs.io). This repository is inspired from [Snakemake Wrapper Repository](https://github.com/snakemake/snakemake-wrappers/). All the Sequana wrappers can can be used in your Snakemake pipelines as well.


[![Tests](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml/badge.svg)](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml)
[![Tests](http://joss.theoj.org/papers/10.21105/joss.00352/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00352)


|||
| --- | --- |
| Overview | A set of Snakemake wappers to help building Sequana pipelines |
|Status | Production |
|Issues | Please fill a report on [github/sequana/sequana-wrappers](https://github.com/sequana/sequana/issues) |
|Python version | Python 3.8, 3.9, 3.10|
|Citation| Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16) ), 352,  [JOSS DOI doi:10.21105/joss.00352 ](http://www.doi2bib.org/bib/10.21105%2Fjoss.00352) |


# Usage

    snakemake --wrapper-prefix https://github.com/sequana/sequana-wrappers

or with a local copy in /home/user:

    cd
    git clone git@github.com:sequana/sequana-wrappers.git sequana_wrappers
    snakemake --wrapper-prefix git+file:///home/user/sequana_wrappers

Note that if you define a variable SEQUANA_WRAPPERS=git+file:///home/user/sequana_wrappers all pipelines will
automatically set the --wraper-prefix to this variable content.

# Notes for developers

## Overview

The wrappers/ directory contains the wrappers themselves. Each sub-directory is dedicated to
a wrapper that is related to a given software/application. A sub directory may have several wrappers (e.g., bwa has a sub directory related to the indexing, and a sub directory related to mapping).

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

### adding a new wrapper in practice

In ./wrappers, add a new wrapper. Copy the existing fastqc wrapper for instance.
Edit the wrapper.py and design a test/Snakefile example for testing. Since you are
a developer, you are problaby developping in a dedicated branch. Let us call it **dev**.

In the test/Snakefile, you should switch from the **main** to the **dev** in the wrapper path:

    wrapper:
        "dev/wrappers/my_new_wrapper"

In order to test your Snakefile, you first need to commit the wrapper.py. Then, execute the Snakefile:

    snakemake -s Snakefile  -j 1 --wrapper-prefix git+file:///YOURPATH/sequana-wrappers/ -f -p

If it fails, edit and commit your wrapper.py and execute again until your Snakefile and wrappers are functional.

Once done, switch back the wrapper path to the **main** branch:

    wrapper:
        "main/wrappers/my_new_wrapper"

Time to include the new wrapper in the continous integration. Go to the root of sequana-wrappers and add a functional test to the end of test.py. Then, test it:

    pytest test.py -k my_new_wrapper -v

You are ready to push and create a pull-requests



## Tagging

You may consider to add a tag. Our convention is to use a tag with the YEAR.MONTH.DAY
where day and month do not include extra zeros. So you would have e.g.::

    23.11.11
    23.2.2

but not 23.02.02



