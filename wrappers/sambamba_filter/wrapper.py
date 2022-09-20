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
from snakemake.shell import shell

# Get rule information (input/output/params...)

input_bam = snakemake.input
output_bam = snakemake.output
options = snakemake.params.get("options", "")
threshold = snakemake.params.get("threshold")
logs = snakemake.log_fmt_shell(stdout=True, stderr=True)

if not threshold:
    raise AssertionError("A threshold must be set in rule params")

shell(
    "sambamba view {options}"
    " --format=bam"
    " --filter='mapping_quality >= {threshold}'"
    " -o {output_bam}"
    " {input_bam} {logs}"
)
