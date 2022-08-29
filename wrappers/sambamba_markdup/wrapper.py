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
options = snakemake.params.get("options", "")
remove_duplicates = snakemake.params.get("remove_duplicates", "False")
tmpdir = snakemake.params.get("tmp_directory")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


if remove_duplicates:
    options += " --remove-duplicates"

if tmpdir:
    options += f" --tmpdir={tmpdir}"

shell("sambamba markdup {input_bam} {output_bam} {options} {log}")
