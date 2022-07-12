#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Dev Team (https://sequana.readthedocs.io)
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
if "alignment" in snakemake.input and "assembly" in snakemake.input:
    assembly = snakemake.input.assembly
    alignment = snakemake.input.alignment
else:
    # alphabetical order
    alignment = snakemake.input[0]
    assembly = snakemake.input[1]


# There is only one output
result = snakemake.output[0]

# default params set to ""
options = snakemake.params.get("options", "")


# the result (FastA sequence) is redirected to stdout by polypolish
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

shell("polypolish {options} {assembly} {alignment} 1>{result} {log}")

