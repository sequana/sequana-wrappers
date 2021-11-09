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

import os

from snakemake.shell import shell

# Get rule information (input/output/params...)

input_bam = snakemake.input[0]
output_bam = snakemake.output[0]
params = snakemake.params
log = snakemake.log



cmd = """sambamba view  {params.options} --format=bam --filter="mapping_quality >= {params.threshold}"       -o {output_bam} {input_bam} 1>{log.out} 2>{log.err}"""

shell(cmd)
