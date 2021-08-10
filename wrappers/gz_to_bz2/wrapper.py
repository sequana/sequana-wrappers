__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import sys
import time
from os import path

from snakemake.shell import shell




# Get directory name
input_file = snakemake.output[0]
output_file= snakemake.output[0]


# check integrity
cmd = "pigz -p{threads} --test {input_file}"
shell(cmd)

# conversion
cmd = "pigz -d -c -p {threads} {input_file} | pbzip2 -p{threads} > {output_file}"
shell(cmd)

# integrity output
cmd = "pbzip2 {output_file} -p{threads} --test"
shell(cmd)

# remove original file
cmd = "rm -f {input_file}"
shell(cmd)

