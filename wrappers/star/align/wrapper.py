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
import math
from pathlib import Path

from snakemake.shell import shell

from sequana import FastA, FastQ



# Get input information to run STAR
input_fastq = snakemake.input.fastq
input_reference = snakemake.input.reference

# Get index from a previously run star_index call
# This serves 2 purposes. First, as a requirment
# for star_align to start, second as an information
# to get the genome directory (prefix) to be used later
input_index = snakemake.input.index
genomeDir = str(Path(input_index).parent)

#
log = snakemake.log_fmt_shell(stdout=True, stderr=True, append=True)


# just an alias
params = snakemake.params

# compulsary parameter
# the BAM file create always finishes with this suffix, meaning the
# first part is the prefix
output = snakemake.output.get("bam")
suffix = "_Aligned.sortedByCoord.out.bam" 
if suffix not in output:
    raise ValueError(f"The output bam file must end with {suffix}")
prefix = output.replace(suffix, "")


# optional documented arguments
readFilesCommand = params.get("readFilesCommand", "zcat")
sjdbOverhang = params.get("sjdbOverhang", 100)
outSAMtype = params.get("outSAMtype", "BAM SortedByCoordinate")
options_first_pass = params.get("options_first_pass", "")
options_second_pass = params.get("options_second_pass", "")
options_indexing = params.get("options_indexing", "")
options = params.get("options", "")

# optional parameters undocumented on purpose
prefix_first_pass = params.get("prefix_first_pass",  f"{prefix}_init_")
prefix_second_pass = params.get("prefix_second_pass", f"{prefix}_")  # note the trailing _
splice_file = params.get("splice_file", f"{prefix}_init_SJ.out.tab")
second_pass_genome_dir = params.get("genome_dir", f"{prefix}_star_2nd_pass")


# Get the optimal number of bases according to STAR documentation
# same as in star_index
ff = FastA(input_reference)
genome_length = sum(ff.lengths)
Nbases = min(14, math.floor(math.log2(genome_length)/2 - 1))
Nbases = params.get("genomeSAindexNbases", Nbases)


# limitGenomeGenerateRAM defaults is 31 000 000 000
# limitBAMsortRAM if 0 will be set to genome index size

# Legacy mode 
if params.get('legacy', False):
    shell(""" echo "Running rna-star 1st pass"
       STAR --genomeDir {genomeDir} \
            --readFilesIn {input_fastq}  \
            --runThreadN {snakemake.threads} \
            --genomeLoad NoSharedMemory \
            --outSAMtype {outSAMtype} \
            --readFilesCommand {readFilesCommand} \
            --outFileNamePrefix {prefix_first_pass} \
           {options_first_pass}  {log}""")

    shell(""" echo "running rna-star genome indexing"
        STAR --runMode genomeGenerate \
             --genomeDir {second_pass_genome_dir} \
             --genomeFastaFiles {input_reference} \
             --sjdbFileChrStartEnd {splice_file} \
             --sjdbOverhang {sjdbOverhang} \
             --runThreadN {snakemake.threads} \
             --genomeSAindexNbases {Nbases} \
             {options_indexing}  {log}""")


    shell(""" echo "Running rna-star 2nd pass"
        STAR --genomeDir {second_pass_genome_dir} \
             --readFilesIn {input_fastq}  \
             --runThreadN {snakemake.threads} \
             --genomeLoad NoSharedMemory \
             --sjdbFileChrStartEnd {splice_file} \
             --outSAMtype {outSAMtype} \
             --readFilesCommand {readFilesCommand} \
             --outFileNamePrefix {prefix_second_pass} \
             {options_second_pass}   {log}""")

# new 2-pass mode recommended by STAR
else:
    # requires 'annotation' field in the input
    cmd = """ echo "Running rna-star 2-pass mode"
        STAR --genomeDir {genomeDir} \
             --twopassMode Basic \
             --twopass1readsN -1\
             --readFilesIn {input_fastq}  \
             --runThreadN {snakemake.threads} \
             --genomeLoad NoSharedMemory \
             --outSAMtype {outSAMtype} \
             --readFilesCommand {readFilesCommand} \
             --outFileNamePrefix {prefix_second_pass} """


    input_annotation = snakemake.input.get("annotation")

    cmd += ' --sjdbGTFfile {input_annotation} --sjdbOverhang {sjdbOverhang} '

    if input_annotation.endswith('gtf'):
        exonParentTranscript = params.get("exonParentTranscript", "transcript_id") 
    elif input_annotation.endswith(('gff', 'gff3')):
        exonParentTranscript = snakemake.params.get("exonParentTranscript", "Parent")
        cmd += ' --sjdbGTFtagExonParentTranscript {exonParentTranscript}'


    cmd += "{options}   {log}"
    shell(cmd)


