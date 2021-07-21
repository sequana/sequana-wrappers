__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"

from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fastq = snakemake.input
params = snakemake.params
log = snakemake.log[0]

# if the content of the file is empty, this will fail. We need to
# touch  a file in such case. 
newinput = []
if isinstance(input_fastq, str):
    input_fastq = [input_fastq]
else:
    input_fastq = input_fastq

if len(input_fastq) != 0:
    for this in input_fastq:
        if this.endswith(".bam") or this.endswith("sam"):
            shell(
            " fastqc -t {snakemake.threads} --outdir {params.wkdir} "
                " {input_fastq} {params.kargs} &> {log}")
        else:
            shell(
                " fastqc -t {snakemake.threads} --outdir {params.wkdir} -f fastq "
                " {input_fastq} {params.kargs} &> {log}")











