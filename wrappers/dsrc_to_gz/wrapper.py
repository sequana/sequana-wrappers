#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2021 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Documentation: http://sequana.readthedocs.io
#  Contributors:  https://github.com/sequana/sequana/graphs/contributors
##############################################################################
import sys
import time
from os import path

from snakemake.shell import shell

# Get directory name
input_file = snakemake.input[0]
output_file = snakemake.output[0]
options = snakemake.params.get("options", "")

cmd = "dsrc d -s -t{snakemake.threads} {options} {input_file} | pigz -p {snakemake.threads} > {output_file}"
shell(cmd)

# Check integrity
cmd = "pigz -p {snakemake.threads} --test {output_file}"
shell(cmd)

# Delete the input file
cmd = "rm -f {input_file}"
shell(cmd)
