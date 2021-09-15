__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

import os
import shutil

from snakemake.shell import shell
from easydev import touch

from sequana_pipetools.snaketools import FileFactory


# Get rule information (input/output/params...)
input_fastq = snakemake.input[0]

# params
params = snakemake.params

from sequana import FastQC, sequana_data
import pylab
pylab.ioff()

ff = FileFactory(input_fastq)

for i, filename in enumerate(ff.realpaths):
    # The ouput files
    output_gc = snakemake.output['gc']
    output_boxplot = snakemake.outoput['boxplot']
    output_json = snakemake.output["json"]

    fastq = FastQC(filename, max_sample=params.max_reads)
    if len(fastq.fastq) != 0:
        pylab.clf()
        fastq.boxplot_quality()
        pylab.savefig(output_boxplot)

        pylab.clf()
        fastq.histogram_gc_content()
        pylab.savefig(output_gc)

        stats = fastq.get_stats()
        stats.to_json(output_json)
    else:
        location = sequana_data("no_data.jpg", "images")
        shutil.copy(location, output_gc)
        shutil.copy(location, output_boxplot)
        # this will be handled inside report_fastq_stats
        touch(output_json)

