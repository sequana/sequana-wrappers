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
assembly = snakemake.input.get("assembly")
alignments = snakemake.input.get("alignments")

if isinstance(alignments, (list, tuple)):
    alignments = " ".join(alignments)

# There is only one output
result = snakemake.output[0]

# default params set to ""
options = snakemake.params.get("options", "")


# the result (FastA sequence) is redirected to stdout by polypolish
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

cmd = "polypolish {options} {assembly} {alignments} 1>{result} {log}"
shell(cmd)

