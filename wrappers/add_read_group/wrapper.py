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

input_bam = snakemake.input[0]
output_bam = snakemake.output[0]
options = snakemake.params.get("options", "")


# we can not have the same options twice
# The one in options (if found) should be kept.
if "-PL" not in options.split():
    PL = snakemake.params.get("PL", "Illumina")
    options += f" -PL {PL}"

if "-LB" not in options.split():
    LB = snakemake.params.get("LB", "unknown")
    options += f" -LB {LB}"

if "-PU" not in options.split():
    PU = snakemake.params.get("PU", "unknown")
    options += f" -PU {PU}"

if "-SM" not in options.split():
    SM = input_bam.strip().rsplit(".", 1)[0].replace(".", "_")
    SM = snakemake.params.get("SM", SM)
    options += f" -SM {SM}"

if "-ID" not in options.split():
    import uuid

    ID = snakemake.params.get("ID", int(uuid.uuid1()))
    options += f" -ID {ID}"


cmd = "picard AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -I {input_bam} -O {output_bam} {options} {log} && samtools index {output_bam}"
shell(cmd)
