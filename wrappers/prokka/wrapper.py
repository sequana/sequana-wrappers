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


assembly = snakemake.input[0]
outfile = Path(snakemake.output[0])
options = snakemake.params.get("options", "")
logs = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    "prokka --force {options}"
    " --cpus {snakemake.threads}"
    " --outdir {outfile.parent}"
    " --prefix {outfile.stem}"
    " {assembly} {logs}"
)