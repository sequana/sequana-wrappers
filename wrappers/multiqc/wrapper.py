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
"""Snakemake wrapper for multiqc."""

__author__ = "Etienne Kornobis, Thomas Cokelaer"
__copyright__ = "Copyright 2021 Sequana Dev team"
__license__ = "BSD"


from os import path

from snakemake.shell import shell

output_dir = path.dirname(snakemake.output[0])
output_name = path.basename(snakemake.output[0])
log = snakemake.log_fmt_shell(stdout=True, stderr=True)


options = snakemake.params.get("options", "")
input_directory = snakemake.params.get("input_directory", ".")
modules = snakemake.params.get("modules", "")
config_file = snakemake.params.get("config_file", "")


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
