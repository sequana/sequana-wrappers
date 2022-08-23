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
import os
import math
from pathlib import Path

from snakemake.shell import shell


# Get input information to run salmon
input_fastq = snakemake.input.fastq
input_index = snakemake.input.index

#
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)

# any salmon options to be provided
options = snakemake.params.get("options", "")

# compulsary parameter. the quant file created by salmon is  always called quant.sf
# actually salmon requires the output directory only
output = snakemake.output[0]
outdir = str(Path(output).parent)

# we also need the directory where is stored the index
indexdir = str(Path(input_index).parent)


# note that we put {optinos} before the -1/-2 and -r options
# because options such as -l A, if provided, must be before 
# the -1/-2 /-r options
if len(input_fastq) == 2:
    # in general, we store results in {sample}/rule_name
    # and results can be e.g. {sample}/rule_name/{sample}.txt
    # but multiqc report for salmon identify in the sample name in the
    # directory where are found the reports. So we must use
    # {sample}/rule_name_{sample}
    shell("""salmon quant -i {indexdir}  {options} \
-1 {input_fastq[0]} -2 {input_fastq[1]} -p {snakemake.threads} \
--validateMappings -o {outdir} {log}""")
else:
    shell("""salmon quant -i {indexdir} {options} \
-r {input_fastq} -p {snakemake.threads} \
--validateMappings -o {outdir}  {log}""")

# in case user gives same name as salmon, need -f to force overwiting
if Path(output).name != 'quant.sf':
    shell("""mv -f {outdir}/quant.sf {output}""")

