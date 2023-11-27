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


input_bam = snakemake.input.bam
input_ref = snakemake.input.ref

output_vcf = snakemake.output[0]

ploidy = snakemake.params["ploidy"]
options = snakemake.params["options"]

log = snakemake.log[0]


cmd = "freebayes {options} --ploidy {ploidy} -f {input_ref} -v {output_vcf} {input_bam} &>{log}"

shell(cmd)
