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

from snakemake import shell

input_bam = snakemake.input[0]
output_bed = snakemake.output[0]
log = snakemake.log_fmt_shell(stdout=False, stderr=True)

try:
    options = snakemake.params["options"]
except AttributeError:
    options = ""


cmd = "samtools depth -m 20000 -aa {input_bam} {options} > {output_bed} {log}"
shell(cmd)
