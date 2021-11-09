import os

from snakemake.shell import shell

from sequana.bedtools import GenomeCov
from sequana.modules_report.coverage import CoverageModule
from sequana.utils import config

# Get rule information (input/output/params...)

input_bed = snakemake.input.bed
input_gbk = snakemake.input.gbk
input_fasta = snakemake.input.fasta

output_html = snakemake.output[0]
report_dir = output_html.rsplit("/", 1)[0]

circular = snakemake.params.get("circular", True)
chunksize = snakemake.params.get("chunksize", 0.5)
double_threshold = snakemake.params.get("double_threshold", 0.5)
gc_size = snakemake.params.get("gc_window_size", 201)
params_k = snakemake.params.get("mixture_models", 2)
high = snakemake.params.get("high_threshold", 4)
low = snakemake.params.get("low_threshold", -4)
window_size = snakemake.params.get("window_size", 3001)



# Run sequana coverage
# read the data
bed = GenomeCov(input_bed, input_gbk, low, high, double_threshold, double_threshold)

# compute GC
bed.compute_gc_content(input_fasta, gc_size, circular)

# we expect the sample name to be the {sample}/something/...
sample = input_bed.split(os.sep)[0]

for chrom in bed:
    if window_size >  len(chrom.df) / 4:
        W = int(len(chrom.df) / 4)
    else:
        W = window_size

    results = chrom.run(W, circular=circular,
                k=params_k, binning=-1, cnv_delta=-1)

    prefix = f"{report_dir}/{chrom.chrom_name}/"
    output_json = prefix + "sequana_summary_coverage.json"

    os.makedirs(prefix, exist_ok=True)
    output_roi = prefix + "rois.csv"
    ROIs = results.get_rois()
    ROIs.df.to_csv(output_roi)
    chrom.plot_coverage(prefix + "coverage.png")

    summary = results.get_summary(caller="sequana_pipeline")
    summary.to_json(output_json)

# Create HTML reports
config.output_dir = report_dir
config.sample_name = sample

CoverageModule(bed)

