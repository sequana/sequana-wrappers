__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import os

from snakemake.shell import shell
from sequana_pipetools import SequanaConfig
from sequana_pipetools.snaketools import DOTParser



# This rule calls snakemake. This is in conflict with the main snakemake call itself.
# Solution: create a new directory where to launch this snakemake

# First, we tried with in a temporary directory but this was creating errs
# most probably because temp dir was handled in the code rather than 
# by snakemake itself. 

# Second, we used a os.chcwd(). Although functional locally, this 
# messes up the main snakemake snakejobs, that could be copied in the 
# new working directory  and then not seen by the main snakemake call
# (latency in the creation of the output files maybe).

# third solution (this one) is to call *cd* the shell commands

input_filename = snakemake.input[0]
output_svg = snakemake.output['svg']
params = snakemake.params

# change relative path to absolute path
def parse_path(dico):
    for key, value in dico.items():
        try:
            if os.path.exists(value):
                dico[key] = os.path.realpath(value)
        # check overflowerror if value is a large int
        except (TypeError, OverflowError):
            try:
                parse_path(value)
            except AttributeError:
                pass
cfg = SequanaConfig(params.configname)
parse_path(cfg.config)
cfg._update_yaml()

cwd = os.getcwd() # if it fails, we must reset the current working directory
try:
    try:
        os.mkdir("rulegraph")
    except:
        pass
    cfg.copy_requirements(target="rulegraph")
    cfg.save(filename="rulegraph" + os.sep + params.configname)
    shell('cd rulegraph; snakemake -s "{input_filename}" --rulegraph --nolock  > rg.dot; cd .. ')
except Exception as err:
    print(err)
    #make sure we come back to the correct directory
    os.chdir(cwd)

# Annotate the dag with URLs
d = DOTParser(cwd + os.sep + "rulegraph/rg.dot")
d.add_urls(mapper=params.mapper)

# Now, create the SVG. Somehow if called dag.svg, this is a conflict
# hence the || true
shell("dot -Tsvg rg.ann.dot  -o {output_svg} || true")
shell("rm -f rg.ann.dot")
