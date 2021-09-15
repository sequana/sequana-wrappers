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
input_fastq = snakemake.fastq

# figure out the expected output directory
done = snakemake.output[0]

# params
params = snakemake.params

from sequana import FastQC, sequana_data
import pylab
pylab.ioff()

ff = FileFactory(input_fastq)

for i, filename in enumerate(ff.realpaths):
    # The ouput files
    formatter = lambda x: params.wkdir + "/" + x.replace(".fastq.gz","")
    output_gc = formatter(ff.basenames[i] + "_gc.png")
    output_boxplot = formatter(ff.basenames[i] + "_boxplot.png")
    output_json = formatter(ff.basenames[i] + ".json")

    fastq = FastQC(filename, max_sample=500000)
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

