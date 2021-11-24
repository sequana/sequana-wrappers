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

from snakemake.shell import shell

# Get rule information (input/output/params...)
input_vcf = snakemake.input.vcf
input_ann = snakemake.input.ann

output_csv = snakemake.output.csv
output_vcf = snakemake.output.vcf
output_html = snakemake.output.html

params = snakemake.params
options = snakemake.params.get("options", "")

log = snakemake.log[0]

# if csv required, add the argument csvStats
if output_csv:
    options = f"{options} -csvStats {output_csv}"


from sequana import SnpEff

mydata = SnpEff(input_ann, log=log)
mydata.launch_snpeff(input_vcf, output_vcf, html_output=output_html, options=options)
