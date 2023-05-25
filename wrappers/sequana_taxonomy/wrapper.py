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

from snakemake import shell

outdir = snakemake.output.html.split("/",1)[0]

cmd = "sequana_taxonomy --file1 {snakemake.input[0]} "
if snakemake.params.paired:
    cmd += " --file2 {snakemake.input[1]} "

confidence = snakemake.params.get('confidence', 0)
level = snakemake.params.get('level', "INFO")

if confidence != 0:
    cmd += " --confidence {confidence} "

# threading
cmd += " --thread {snakemake.threads} --databases"

# handle databases
for this in snakemake.params['databases']:
    if this != "toydb":
        assert os.path.exists(this), f"databases {this} does not exits"
    cmd += f" {this} "

cmd += " --output-directory {outdir} "
cmd += " --level {level} "

# we will save the unclassified files whatsover
if snakemake.params.store_unclassified:
    if snakemake.params.paired:
        cmd += " --unclassified-out unclassified#.fastq"  # this syntax replaces # with 1 and 2
    else:
        cmd += " --unclassified-out unclassified.fastq"

# add the option at the end to overwrite previous choices
cmd += snakemake.params.get("options", "")

# even tough we do not save the unclassified, we want to create a dummy 
# file since this ouptut file is expected.
cmd += f" && touch {snakemake.output.unclassified}"
shell(cmd)


