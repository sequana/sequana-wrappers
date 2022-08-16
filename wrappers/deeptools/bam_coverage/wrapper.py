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

bam = snakemake.input[0]
outfile = snakemake.output[0]
options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

bam_index = Path(bam + ".bai")
if not bam_index.exists():
    shell("samtools index {bam}")

shell(
    "bamCoverage --bam {bam}"
    " --outFileName {outfile}"
    " --numberOfProcessors {snakemake.threads}"
    " {options} {log}"
)
