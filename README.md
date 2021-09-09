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


A wrappers contains these files:

- environment.yaml to tell what are the required packages to be installed for the test (and wrapper) to work.
- wrappers.py the wrapper itself. There is no specific instructions here except to write good code as much as possible :-)
- example.sml is used to document your wrapper. It is an example of a Snakefile using the wrapper. It is not intended to be tested for now. 
  Instead it is meant to have a docstring that will be used for documentation.
* a test/ directory with one file called **Snakefile** 

So to add a new wrappers:

1. In the .wrappers/ directory create a directory **in small caps** (e.g. for this example let us call it **example**
2. in ./wrappers/example, add the files mentionned above
3. in the root of the repository, there is a test.py file. Please see the bottom of this file and add a test for your wrapper.
4. Test the new wrapper:

   pytest test.py -k test_example
   
tooltips: if you use a branch named add_wrapper, the Snakefile you have just added should use the wrapper your_branch/wrappers/example temporarily until you push the new wrapper.

## Notes about the config file used in sequana

- Rules may use the threads key, which can be defined in the general config file.
- the **params** section should contain a key called**options** also define in the config file of most pipelines.
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


## Faqs

When adding new recipe / testing, you may face several issues. 

Make sure you have added/commited the files you want to test.

If you get this message::

    Failed to open source file /tmp/tmpw9y6n0ju/main/wrappers/bowtie1/build/environment.yaml


To test locally, create a branch, add the recipe. The wrapper must replace the
'main' by your branch name and then you can test the wrapper locally as
follows::

   cd wrapper/your_recipe/test
   snakemake -s Snakefile  -j 1 --wrapper-prefix git+file:///YOURPATH/sequana-wrappers/ -f -p

Then, once it works, do not forget to replace the branch name in the test.Snakefile with "main"





