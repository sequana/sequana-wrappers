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

input_fasta = snakemake.input[0]
input_ann = snakemake.input[1]

output_fasta = snakemake.output[0]

log = snakemake.log[0]

options = snakemake.params.get("options", "")

# real stuff is here:
from sequana import SnpEff

if input_ann.endswith(".gbk"):
    snpeff = SnpEff(input_ann, log=log, build_options=options)
elif input_ann.endswith("gff") or input_ann.endswith("gff3"):
    snpeff = SnpEff(input_ann, log=log, fastafile=input_fasta, build_options=options)
else:
    raise IOError("Your annotation file does not end with gbk or gff/gff3 extension")

snpeff.add_locus_in_fasta(input_fasta, output_fasta)
