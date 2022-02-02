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

options = snakemake.params.get("options", "")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

n = len(snakemake.input.fastq)

if n == 1:
    reads = "-U {}".format(*snakemake.input.fastq)
elif n == 2:
    reads = "-1 {} -2 {}".format(*snakemake.input.fastq)
else:
    raise ValueError("input->fastq must have 1 (single-end) or 2 (paired-end) elements.")

shell(
    "(bowtie2 --threads {snakemake.threads} {options} "
    "-x {snakemake.params.index} {reads} "
    "| samtools view -Sbh -o {snakemake.output.bam} -) {log}"
)

shell("samtools sort {snakemake.output.bam} > {snakemake.output.sorted}")
shell("samtools index {snakemake.output.sorted}")
