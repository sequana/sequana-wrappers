from snakemake.shell import shell

left = snakemake.input.get("left")
assert left is not None, "input-> left is a required input parameter"
right = snakemake.input.get("right")
if right:
    input_cmd = f" --left {left} --right {right}"
else:
    input_cmd = f" --single {left}"

# Check user input but use kallisto by default
est_method = snakemake.params.get("est_method", "kallisto")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell(
    """align_and_estimate_abundance.pl \
    --seqType fq \
    {input_cmd} \
    --transcripts {snakemake.input.fasta} \
    --est_method {est_method} \
    --trinity_mode \
    --SS_lib_type {snakemake.params.lib_type} \
    --thread_count {snakemake.threads} \
    --output_dir {snakemake.output.outdir} \
    --prep_reference {log} \
    """
)
