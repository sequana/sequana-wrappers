#
#  This file is part of Sequana software
#
#  Copyright (c) 2016-2022 - Sequana Dev Team (https://sequana.readthedocs.io)
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

# Get rule information (input/output/params...)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# alias to input files
inputs = snakemake.input.inputs

# outputs need some tuning, let us make sure we are handling a list (could be a filename or list of filenames)
controls = snakemake.input.controls
#controls = [controls] if isinstance(controls, str) else controls

# same for the outputs
outputs = snakemake.output
outputs = [outputs] if isinstance(outputs, str) else outputs


# all output filenames must have the same prefix !
prefixes = [x if Path(x).is_dir() else Path(x).parent for x in outputs]
if len(prefixes) != 1:
    prefixes = ",".join(prefixes)
    raise ValueError(f"Sequana-wrappers:macs3 error. output files must all have the same prefixes. Found {prefixes}")

# from the output files, let us figure out the output prefixes.
if Path(outputs[0]).is_dir():
    # if no  / separator provided, then outdir is local . directory
    outdir = outputs[0]
else:
    # we remove the last file, the remaining prefix is the outdir directory
    outdir = Path(outputs[0]).parent

# ---------------------------------------- manage the input parameters and create the command
params = snakemake.params

qvalue = params.get("qvalue", 0.05)
options = params.get("options", " --keep-dup all ")

# -B --SPMR is to save fragment pileup useful to save into bigwig later on
cmd = (
    "macs3 callpeak -B --SPMR "
    " -t {inputs} "
    " -c {controls} "
    " -g {params.genome_size} "
    " -n {params.prefix} "
    "--bw {params.bandwidth} "
    " {options} -q {qvalue} "
)

# is it paired data or not ?
if params.paired:
    cmd += " -f BAMPE "
else:
    cmd += " -f BAM "

# switch between the narrow and broad cases
if params.mode == "narrow":
    shell(cmd + " --outdir {outdir} {log}")
elif params.mode == "broad":
    shell(cmd + " --outdir {outdir} --broad --broad-cutoff {params.broad_cutoff} {log}")
else:
    raise ValueError(f"Sequana-wrappers:macs3 error. The 'mode' key in the params section can be only set to 'narrow' or 'broad'; you gave {params.mode}")

