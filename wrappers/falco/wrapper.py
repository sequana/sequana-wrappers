__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import os

from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fastq = snakemake.input
params = snakemake.params
log = snakemake.log[0]

# if the content of the file is empty, this will fail. We need to
# touch  a file in such case.
from sequana import FastQ

# Not that falco v0.1 process all fastq sequentially and erase the previous one
# with new ones. So we process the first one only. 

fastq = FastQ(input_fastq[0])
if len(fastq) == 0:
    pass
else:
   shell("""falco -t {snakemake.threads} --outdir {params.working_directory} {input_fastq[0]} {params.options} &> {log}""")




