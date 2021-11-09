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
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

ID = snakemake.params.get(ID)
LB = snakemake.params.get(LB, "unknown")
PL = snakemake.params.get(PL, "ILLUMINA")
PU = snakemake.params.get(PU, "unknown")
SM = snakemake.params.get(SM)

RG = f"{ID={ID} LB={LB} PL={PL} PU={SM} SM={SM}"


cmd = "picard AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT
I={snakemake.input[0]} O={snakemake.output.bam} {RG}  {log} && samtools index {snakemake.output.bam}"

shell(cmd)

