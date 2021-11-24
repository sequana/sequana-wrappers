"""Snakemake wrapper for BUSCO assessment"""

__author__ = "Tessa Pierce"
__copyright__ = "Copyright 2018, Tessa Pierce"
__email__ = "ntpierce@gmail.com"
__license__ = "MIT"

from snakemake.shell import shell
from os import path


log = snakemake.log_fmt_shell(stdout=True, stderr=True)
options = snakemake.params.get("options", "")
mode = snakemake.params.get("mode")
assert mode is not None, "please input a run mode: genome, transcriptome or proteins"
lineage = snakemake.params.get("lineage")
assert lineage is not None, "please input the path to a lineage for busco assessment"

stripped_output = snakemake.output[0].rstrip("/")
out = path.basename(stripped_output)
out_dirname = path.dirname(stripped_output)
out_path = " --out_path {} ".format(out_dirname) if out_dirname else ""
# change short summary filename for convenience and to have data split in multiqc report
try:
    move_file = (
        f"&& mv {out_dirname}/{out}/short_summary.*.txt {out_dirname}/{out}/{snakemake.params.short_summary_filename}"
    )
except AttributeError:
    move_file = ""

download_path_dir = snakemake.params.get("download_path", "")
download_path = " --download_path {} ".format(download_path_dir) if download_path_dir else ""

# note: --force allows snakemake to handle rewriting files as necessary
# without needing to specify *all* busco outputs as snakemake outputs
shell(
    "busco --in {snakemake.input} --out {out} --force "
    "{out_path} "
    "--cpu {snakemake.threads} --mode {mode} --lineage {lineage} "
    "{download_path} "
    "{options} {log} {move_file}"
)
