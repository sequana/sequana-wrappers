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
from pathlib import Path
from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fasta = snakemake.input.get("fasta", snakemake.input[0])
input_gff = snakemake.input.get("gff", "")

# output use to get working dir. There is just one output, so
# let us use [0]
done = snakemake.output[0]

# retrieve a default wkdir from the output
wkdir = str(Path(done).parent)

# if user provide a wkdir, overwrite it
wkdir = snakemake.params.get("wkdir", wkdir)

# options
options = snakemake.params.get("options", "")

# log file
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# actual code
shell(f"gffread {input_gff} -g {input_fasta} -w {wkdir}/salmon_transcript.fa {log} ")
shell(f"salmon index -p {snakemake.threads} -t {wkdir}/salmon_transcript.fa -i {wkdir} {options} {log} ")

# trigger file
shell("touch {done}")

