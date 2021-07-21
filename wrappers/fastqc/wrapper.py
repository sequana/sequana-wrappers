__author__ = "Thomas Cokelaer"
__copyright__ = "Copyright 2021, Sequana Team"
__email__ = "thomas.cokelaer@pasteur.fr"
__license__ = "BSD-3"


from snakemake.shell import shell


# Get rule information (input/output/params...)
input_fastq = snakemake.input

# figure out the expected output directory
done = snakemake.input[0]
outdir = done.rsplit("/")[0]

params = snakemake.params
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


# Note that if the input file is empty, fastqc will fail. 
# We need to touch  a file in such case and your Snakefile
# should define as output a "done" file.

if isinstance(input_fastq, str):
    input_fastq = [input_fastq]
else:
    input_fastq = input_fastq

if len(input_fastq) != 0:
    for fastq_file in input_fastq:
        if fastq_file.endswith((".bam", "sam")):
            shell(
            " fastqc -t {snakemake.threads} --outdir {outdir} "
                " {fast_file} {params.options} &>> {log}")
        else:
            shell(
                " fastqc -t {snakemake.threads} --outdir {outdir} -f fastq "
                " {fastq_file} {params.options} &>> {log}")

shell("touch({done})")

