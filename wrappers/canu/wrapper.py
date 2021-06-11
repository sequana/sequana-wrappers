__author__ = "Dimitri Desvillechabrol"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "dimitri.desvillechabrol@pasteur.fr"
__license__ = "BSD-3"

import sys
import time
from os import path

from snakemake.shell import shell


suffix_dict = {
    '-correct': '.correctedReads.fasta.gz',
    '-trim': '.trimmedReads.fasta.gz',
}

options = snakemake.params.get('options', '')
use_grid = snakemake.params.get('use_grid', False)
step = snakemake.params.get('step', '')

# Get directory name
output_file = snakemake.output[0]
output_dir = path.dirname(output_file)
# Get prefix name
filename = path.basename(output_file)
prefix = filename.replace(suffix_dict.get(step, '.contigs.fasta'), "")

# remove previous done/failed files
shell(f"rm -f {output_dir}/canu{step}.done {output_dir}/canu{step}.failed")

# run canu
shell(
    f"canu {step}"
    f" -p {prefix}"
    f" -d {output_dir}"
    " genomeSize={snakemake.params.genome_size}"
    " maxThreads={snakemake.threads}"
    f" useGrid={use_grid}"
    f" {options}"
    f" onSuccess='touch canu{step}.done'"
    f" onFailure='touch canu{step}.failed'"
    " -pacbio {snakemake.input[0]}"
)

# wait canu if grid is used
if use_grid:
    while not path.exists(f'{output_dir}/canu{step}.done'):
        time.sleep(60)
        if path.exists(f'{output_dir}/canu{step}.failed'):
            sys.exit(1)            
