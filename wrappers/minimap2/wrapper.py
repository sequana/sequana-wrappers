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
input_fastq = snakemake.input["fastq"]
input_reference = snakemake.input["reference"]
output = snakemake.output
params = snakemake.params
threads = snakemake.threads


shell(
    """
minimap2 -t {threads} {input_reference} {input_fastq} {params.options} -a | samtools sort -@ {threads} -o {output}
"""
)
