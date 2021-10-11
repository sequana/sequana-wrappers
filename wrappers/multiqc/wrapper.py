"""Snakemake wrapper for multiqc."""

__author__ = "Etienne Kornobis, Thomas Cokelaer"
__copyright__ = "Copyright 2021 Sequana Dev team"
__license__ = "BSD"


from os import path

from snakemake.shell import shell

print("--------------------")
output_dir = path.dirname(snakemake.output[0])
output_name = path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


options = snakemake.params.options
input_directory = snakemake.params.input_directory
modules = snakemake.params.modules
config_file = snakemake.params.config_file


print(options)
print(input_directory)
# if config file not provided, should be set to empty string
if config_file.strip():
    config_file = f" -c {config_file}"

# if modules is not provided, should be set to empty string

module_options = ""
if modules.strip():
    modules = modules.split()
    for module in modules:
        module_options += f" -m {module}"


shell(
    "multiqc"
    " {input_directory}"
    " --force"
    " {options}"
    " {module_options}"
    " -o {output_dir}"
    " -n {output_name}"
    " {config_file}"
    " {log}"
)
