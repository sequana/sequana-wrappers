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
from pathlib import Path

from snakemake.shell import shell


fastq = snakemake.input
contigs = snakemake.output
logs = snakemake.log_fmt_shell(stdout=True, stderr=True)

is_paired = len(fastq) == 2
input_file = f"-1 {fastq[0]} -2 {fastq[1]}" if is_paired else f"-s {fastq[0]}"

options = snakemake.params.get("options", "")
mode = snakemake.params.get("mode", "normal")
outdir = Path(str(contigs)).parent

shell(
    "unicycler --mode {mode}"
    " --threads {snakemake.threads}"
    " {input_file}"
    " -o {outdir}"
    " {options} {logs}"
    " && cp {outdir}/assembly.fasta {contigs}"
)
