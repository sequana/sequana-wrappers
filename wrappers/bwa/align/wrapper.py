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
__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Dev Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD"


from snakemake.shell import shell
from shutil import which 

fastq = snakemake.input.fastq
reference = snakemake.input.reference
options = snakemake.params.get("options", "")
tmpdir = snakemake.params.get("tmp_directory", "")
output_sorted_bam = snakemake.output.sorted
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


if which("pbwa"):
    bwa_exe = "pbwa"
else:
    bwa_exe = "bwa"

sambamba_tmp = f"--tmpdir={tmpdir}" if tmpdir else ""


shell(
    "({bwa_exe} mem -t {snakemake.threads} {options} {reference} {fastq}"
    " | sambamba view -t {snakemake.threads} -S -f bam -o /dev/stdout /dev/stdin"
    " | sambamba sort /dev/stdin -o {output_sorted_bam} -t {snakemake.threads} {sambamba_tmp}"
    ") {log}"
)
