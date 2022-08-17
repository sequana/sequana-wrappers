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

from sequana import FastA
import math

# Get rule information (input/output/params...)
input_fasta = snakemake.input.get("fasta", snakemake.input[0])
input_annotation = snakemake.input.get("annotation", "")

# The options. we specify sjdbOverhang manually since it should
# be important if annotation are provided.
options = snakemake.params.get("options", "")
sjdbOverhang = snakemake.params.get("sjdbOverhang", 100)
wkdir = snakemake.params.get("wkdir", 'star_index')

# log file
log = snakemake.log_fmt_shell(stdout=True, stderr=True)



# Get the optimal number of bases according to STAR documentation:
ff = FastA(input_fasta)
genome_length = sum(ff.lengths)
Nbases = min(14, math.floor(math.log2(genome_length)/2 - 1))
Nbases = snakemake.params.get("genomeSAindexNbases", Nbases)


cmd = f"STAR --runMode genomeGenerate --genomeFastaFiles {input_fasta} --genomeDir {wkdir} --runThreadN {snakemake.threads} --genomeSAindexNbases {Nbases} "

if input_annotation:

    # compulsary input file
    cmd += ' --sjdbGTFfile {input_annotation} --sjdbOverhang {sjdbOverhang} '
    if input_annotation.endswith('gtf'):
        # default star option for GTF
        exonParentTranscript = snakemake.params.get("exonParentTranscript", "transcript_id")
    elif input_annotation.endswith(('gff', 'gff3')):
        # default star option for GFF
        exonParentTranscript = snakemake.params.get("exonParentTranscript", "Parent")
        cmd += ' --sjdbGTFtagExonParentTranscript {exonParentTranscript}'


cmd += " {options} {log} "

shell(cmd)

# figure out the expected output file and create/touch it 
done = snakemake.output[0]
shell("touch {done}")

