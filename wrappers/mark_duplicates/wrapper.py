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
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

input_bam = snakemake.input[0]
output_bam = snakemake.output[0]
output_metrics = snakemake.output[1]

options = snakemake.params.get("options", "")
params = snakemake.params


cmd = (
    "picard MarkDuplicates  I={input_bam} O={output_bam} M={output_metrics} "
    " REMOVE_DUPLICATES={params.remove_dup} TMP_DIR={params.tmpdir}  && "
    " samtools index {output_bam} {options} {log} && samtools index {output_bam}"
)
shell(cmd)


