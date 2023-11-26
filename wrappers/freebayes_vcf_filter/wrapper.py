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

from sequana.modules_report.coverage import CoverageModule
from sequana.utils import config

# Get rule information (input/output/params...)


input_vcf = snakemake.input[0]

output_vcf = snakemake.output.vcf
output_csv = snakemake.output.csv
output_html = snakemake.output.html

params = snakemake.params

# use the sample name based on output HTML file otherwise, user can provide one
report_dir = params.get("report_dir", output_html.rsplit("/", 1)[0])


from sequana.freebayes_vcf_filter import VCF_freebayes
from sequana.modules_report.variant_calling import VariantCallingModule
from sequana.utils import config

v = VCF_freebayes(input_vcf)
filter_v = v.filter_vcf(params["filter_dict"])
filter_v.to_vcf(output_vcf)
filter_v.to_csv(output_csv)

# the HTML report
config.output_dir = report_dir
VariantCallingModule(filter_v)
