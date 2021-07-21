__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"


from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fastq = snakemake.input

# figure out the expected output directory
done = snakemake.output[0]

# figure out the output directory (just removing the filename.done)
outdir = done.rsplit("/", 1)[0]

params = snakemake.params
# we do not use this alias because we will concatenate possibly 2 output of
# fastqc into a single log file for the paired data
# log = snakemake.log_fmt_shell(stdout=True, stderr=True)
log = snakemake.log[0]


# Note that if the input file is empty, fastqc creates a HTML file

# ( ͠° ͟ʖ ͡°) isinstance
if isinstance(input_fastq, str):
    input_fastq = [input_fastq]
else:
    input_fastq = input_fastq

for fastq_file in input_fastq:
    if fastq_file.endswith((".bam", "sam")):
        shell(
        " fastqc -t {snakemake.threads} --outdir {outdir} "
            " {fast_file} {params.options} &>> {log}")
    else:
        shell(
            " fastqc -t {snakemake.threads} --outdir {outdir} -f fastq "
            " {fastq_file} {params.options} &>> {log}")

shell("touch {done} ")

