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


fastq = snakemake.input.fastq
assembly = snakemake.input.assembly
outfile = Path(snakemake.output[0])
options = snakemake.params.get("options", "")
preset = snakemake.params.get("preset", "single")
logs = snakemake.log_fmt_shell(stdout=True, stderr=True)

reference_options = ""
if "reference" in snakemake.params:
    reference_options += f"-r {snakemake.params.reference}"
    if "annotation" in snakemake.params:
        reference_options += f" --feature {snakemake.params.annotation}"

if len(fastq) == 2:
    input_fastq = f"-1 {fastq[0]} -2 {fastq[1]}"
else:
    input_fastq = f"--{preset} {fastq}"

outdir = outfile.parent

shell(
    "quast.py {options}"
    " -t {snakemake.threads}"
    " {reference_options}"
    " {assembly}"
    " {input_fastq}"
    " --output-dir {outdir} {logs}"
    " && touch {outfile}"
)
