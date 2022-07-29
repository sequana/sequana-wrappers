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


fastq = snakemake.input.get("fastq")
assembly = snakemake.input.get("assembly")
options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

model = snakemake.params.get("model")
if not model:
    raise AssertionError("Please input the path to a model.")

# get output directory
output_file = Path(snakemake.output[0])
output_dir = output_file.parent
default_file = output_dir / "consensus.fasta"

shell(
    "(medaka_consensus {options}" 
    " -t {snakemake.threads}"
    " -i {fastq}"
    " -d {assembly}"
    " -o {output_dir}"
    " -m {model}"
    " && mv {default_file} {output_file}) {log}"
)
