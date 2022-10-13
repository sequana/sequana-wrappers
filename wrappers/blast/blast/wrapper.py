#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Dev Team (https://sequana.readthedocs.io)
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  Website:       https://github.com/sequana/sequana
#  Documentation: http://sequana.readthedocs.io
#  Contributors:  https://github.com/sequana/sequana/graphs/contributors
##############################################################################
from pathlib import Path
from os import path

from snakemake.shell import shell


blastdb = snakemake.input.blastdb
query = snakemake.input.query
blast_type = snakemake.params.blast_type  # required params
db_type = snakemake.params.get("db_type", "")
outfmt = snakemake.params.get("outfmt")
evalue = snakemake.params.get("evalue")

log = snakemake.log_fmt_shell(stdout=False, stderr=True)

try:
    blastdb = Path(blastdb)
    if db_type:
        blastdb = blastdb / db_type
except TypeError:
    # handle multiext case that create a list of input file
    blastdb, _ = path.splitext(blastdb[0])

outfmt = f"-outfmt '{outfmt}'" if outfmt else ""
evalue = f"-evalue {evalue}" if evalue else ""

shell(
    "{blast_type}"
    " {snakemake.params.options}"
    " -query {snakemake.input.query}"
    " {outfmt}"
    " {evalue}"
    " -db {blastdb}"
    " -num_threads {snakemake.threads}"
    " -out {snakemake.output} {log}"
)
