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

cmd = "dsrc d -s -t{threads} {snakemake.params} {input_file} | pigz -p {threads} > {output_file}"
shell(cmd)

# Check integrity
cmd = "pigz -p {threads} --test {output_file}" 
shell(cmd)

# Delete the input file
cmd = "rm -f {input_file}"
shell(cmd)

