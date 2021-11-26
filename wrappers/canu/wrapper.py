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
__author__ = "Dimitri Desvillechabrol"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "dimitri.desvillechabrol@pasteur.fr"
__license__ = "BSD-3"

import sys
import time
from os import path

from snakemake.shell import shell


suffix_dict = {
    "-correct": ".correctedReads.fasta.gz",
    "-trim": ".trimmedReads.fasta.gz",
}

options = snakemake.params.get("options", "")
use_grid = snakemake.params.get("use_grid", False)
step = snakemake.params.get("step", "")
preset = snakemake.params.get("preset", "pacbio")
# canu keeps maxThreads even with use_grid=True ( ͠° ͟ʖ ͡°)
max_threads = "" if use_grid else f"maxThreads={snakemake.threads}"

# Get directory name
output_file = snakemake.output[0]
output_dir = path.dirname(output_file)
# Get prefix name
filename = path.basename(output_file)
prefix = filename.replace(suffix_dict.get(step, ".contigs.fasta"), "")

# remove previous done/failed files
shell(f"rm -f {output_dir}/canu{step}.done {output_dir}/canu{step}.failed")

# run canu
shell(
    "canu {step}"
    " -p {prefix}"
    " -d {output_dir}"
    " genomeSize={snakemake.params.genome_size}"
    " {max_threads}"
    " useGrid={use_grid}"
    " {options}"
    " onSuccess='touch canu{step}.done'"
    " onFailure='touch canu{step}.failed'"
    " -{preset} {snakemake.input[0]}"
)

# wait canu if grid is used
if use_grid:
    while not path.exists(f"{output_dir}/canu{step}.done"):
        time.sleep(60)
        if path.exists(f"{output_dir}/canu{step}.failed"):
            sys.exit(1)
