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
input_fastq = snakemake.input
params = snakemake.params
log = snakemake.log[0]

# if the content of the file is empty, this will fail. We need to
# touch  a file in such case.
from sequana import FastQ

# Not that falco v0.1 process all fastq sequentially and erase the previous one
# with new ones. So we process the first one only.

fastq = FastQ(input_fastq[0])
if len(fastq) == 0:
    pass
else:
    shell(
        """falco -t {snakemake.threads} --outdir {params.working_directory} {input_fastq[0]} {params.options} &> {log}"""
    )
