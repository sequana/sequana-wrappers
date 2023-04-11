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


prefix = snakemake.input.index.replace(".1.ebwt", "")
options = snakemake.params.get("options", "")
threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# fastq could be a single filename or a list
reads = snakemake.input.fastq

try:
    # if a single string is provided, this is a Single-End filename
    # we just use strip to fall back
    reads = reads.strip()
    reads = f" {reads}"
except AttributeError:
    reads = f" -1 {reads[0]} -2 {reads[1]}"


shell(
    "(bowtie -S {options} -p {threads} -x {prefix} {reads}"
    "| samtools view -Sbh  - > {snakemake.output.bam}) {log}"
)


try:
    snakemake.output.sorted
    # sort the bam
    shell("samtools sort -@ {threads} -o {snakemake.output.sorted} {snakemake.output.bam}")
    # and index it
    shell("samtools index {snakemake.output.sorted}")
except AttributeError:
    # FIXME. could add a logger.warning here possibly in the future
    pass

