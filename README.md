# The Sequana Wrapper Repository

The Sequana Wrapper Repository is a collection of reusable wrappers used in Sequana Pipelines inspired from [Snakemake Wrapper Repository](https://github.com/snakemake/snakemake-wrappers/).




[![Tests](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml/badge.svg)](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml)
[![Tests](http://joss.theoj.org/papers/10.21105/joss.00352/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00352)


|||
| --- | --- |
| Overview | A set of tools to help building or using sequana pipelines |
|Status | Production |
|Issues | Please fill a report on [github/sequana/sequana-wrappers](https://github.com/sequana/sequana/issues) |
|Python version | Python 3.7, 3.8|
|Citation| Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16) ), 352,  [JOSS DOI doi:10.21105/joss.00352 ](http://www.doi2bib.org/bib/10.21105%2Fjoss.00352) |


# Notes for developers

The wrappers directory contains the wrappers. Each sub-directory is dedicated to
a software. In the wrapper/software sub-directory you may have several
sub-directories. For example, the wrappers/bowtie1 directory contains two
wrappers called **align** and **build**.

A wrapper directory must contain a file called **wrapper.py** where the
developers must provide the core of the wrapper. There is no specific
instructions here except to write good code as much as possible.

A wrapper directory should have a **test** directory for continuous integration
with a **Snakefile** to be tested. The test must be added in the **test.py**
(test.py)[test.py]. 

For testing purposes, you should also add a file called **environment.yaml** 
to tell what are the required packages to be installed for the test (and wrapper) 
to work.

Finally, for documentation, we kindly ask the developer to create a **README.md** described here below. 

So to add a new wrapper called **example**:

1. Create a directory wrappers/example/ **in small caps** (e.g. for this example let us call it **example**
2. in ./wrappers/example, add the files mentionned above (wrappers.py,
   environment.yaml and README.md).
3. Create a test called **Snakefile** in ./wrappers/example/test with a workable
   example.
4. Add the test code in the file test.py based on previous examples.
4. Test the new wrapper:

   pytest test.py -k test_example

## Notes about the config file used in sequana

- Rules may use the threads key, which can be defined in the general config file.
- the **params** section should contain a key called **options** also define in the config file of most pipelines.
- keys or parameters related to directories and files should use the *_directory* or *_file* suffices

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

## Documentation

Please see the wrappers/fastqc/README.md example. The file must be in markdown
format. It must contain a **Documentation** and **Example** sub sections. If a
**Configuration** section is found, it is also added to the documentation to be
found in https://sequana.readthedocs.io



## Faqs

When adding new recipe / testing, you may face several issues. 

Make sure you have added/commited the files you want to test.

Generally, you should create a branch, add the recipe. In the test file, when setting the wrapper
you should replace the first 'main' by your branch name and then you can test the wrapper locally as
follows::

   cd wrapper/your_recipe/test
   snakemake -s Snakefile  -j 1 --wrapper-prefix git+file:///YOURPATH/sequana-wrappers/ -f -p

Once it works, do not forget to replace the branch name in the test.Snakefile with "main"





