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

input_gff = snakemake.input["gff"]
input_bam = snakemake.input["bam"]
options = snakemake.params.get("options")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# add the strandness
strandness = snakemake.params.get("strandness")
if str(strandness).strip() not in ["0", "1", "2"]:
    raise ValueError("the strandness must be set to 0, 1, or 2")
options += f" -s {strandness} "

if not "-F SAF" in options:
    # add required feature and attribute options if not in SAF mode
    options += f" -t {snakemake.params.feature} "
    options += f" -g {snakemake.params.attribute} "

shell(
    """featureCounts -T {snakemake.threads} {options} -a {input_gff} -o {snakemake.output.counts} {input_bam} {log}"""
)
