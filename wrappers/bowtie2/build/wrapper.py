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

from snakemake.shell import shell

options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
indexbase = snakemake.output[0].replace(".1.bt2", "")

shell("bowtie2-build --threads {snakemake.threads} {options} {snakemake.input.reference} {indexbase} {log}")
shell("samtools faidx {snakemake.input.reference}")
