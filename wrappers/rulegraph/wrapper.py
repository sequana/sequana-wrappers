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
from pathlib import Path

from snakemake.shell import shell
from sequana_pipetools import SequanaConfig
from sequana_pipetools.snaketools import DOTParser


# This rule calls snakemake. This is in conflict with the main snakemake call itself.
# Solution: create a new directory where to launch this snakemake

# First, we tried with in a temporary directory but this was creating errs
# most probably because temp dir was handled in the code rather than
# by snakemake itself.

# Second, we used a os.chcwd(). Although functional locally, this
# messes up the main snakemake snakejobs, that could be copied in the
# new working directory  and then not seen by the main snakemake call
# (latency in the creation of the output files maybe).

# third solution (this one) is to call *cd* the shell commands

input_filename = snakemake.input[0]
output_dot = snakemake.output[0]

# get params
config_filename = snakemake.params.get("configname")
mapper = snakemake.params.get("mapper", {})
required_local_files = snakemake.params.get("required_local_files", [])


# change relative path to absolute path
def parse_path(dico: dict):
    for key, value in dico.items():
        try:
            if value:
                try:
                    dico[key] = str(Path(value).resolve(strict=True))
                except FileNotFoundError:
                    pass
        # check overflowerror if value is a large int
        except (TypeError, OverflowError):
            try:
                parse_path(value)
            except AttributeError:
                pass


def link_required_files(required_paths: list):
    """Symbolic link of required files"""
    for path in required_paths:
        try:
            if Path(path).exists():
                Path(f"rulegraph/{path}").symlink_to(f"../{path}")
            else:
                raise FileNotFoundError(f"Required local file {path} is not found.")
        except FileExistsError:
            pass


cfg = SequanaConfig(config_filename)
parse_path(cfg.config)
cfg._update_yaml()

cwd = Path.cwd()  # if it fails, we must reset the current working directory
try:
    Path("rulegraph").mkdir(exist_ok=True)
    if config_filename:
        cfg.save(filename=Path("rulegraph") / config_filename)
    link_required_files(required_local_files)
    shell('cd rulegraph && snakemake -s "{input_filename}" --rulegraph --nolock  > rg.dot; cd ..')
except Exception as err:
    print(err)
    # make sure we come back to the correct directory
    os.chdir(cwd)

# Annotate the dag with URLs
d = DOTParser(cwd / "rulegraph" / "rg.dot")
d.add_urls(output_filename=output_dot, mapper=mapper)

