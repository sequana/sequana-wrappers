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
from pathlib import Path

from snakemake.shell import shell


fastq = snakemake.input
contigs = snakemake.output.contigs
scaffolds = snakemake.output.scaffolds
logs = snakemake.log_fmt_shell(stdout=True, stderr=True)

is_paired = len(fastq) == 2
input_file = f"-1 {fastq[0]} -2 {fastq[1]}" if is_paired else f"-s {fastq[0]}"


kmers = snakemake.params.get("k", "'auto'")
options = snakemake.params.get("options", "")
preset = snakemake.params.get("preset")
memory = snakemake.params.get("memory", "32")

if preset.startswith("meta") and not is_paired:
    raise NotImplementedError("Cannot use meta SPAdes with single-end data")

preset_option = f"--{preset}" if preset else ""
outdir = Path(scaffolds).parent

shell(
    "spades.py {preset_option}"
    " -k {kmers}"
    " --memory {memory}"
    " --threads {snakemake.threads}"
    " {input_file}"
    " -o {outdir}"
    " {options} {logs}"
    " && cp {outdir}/contigs.fasta {contigs}"
    " && cp {outdir}/scaffolds.fasta {scaffolds}"
)
