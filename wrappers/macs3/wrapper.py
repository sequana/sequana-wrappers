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
from snakemake.shell import shell

# Get rule information (input/output/params...)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# convert IP and control input into list (could be a single filename or list of filenames)
ip_bam = [snakemake.input.IP] if isinstance(snakemake.input.IP, str) else snakemake.input.IP
control_bam = [snakemake.input.control] if isinstance(snakemake.input.control, str) else snakemake.input.control


# what is the output
outputs = [snakemake.output] if isinstance(snakemake.output, str) else snakemake.output


# all output filenames must have the same prefix !
prefixes = [x.rsplit("/",1)[0] for x in outputs]
if len(prefixes) != 1:
    raise ValueError(f"Sequana-wrappers:macs3 error. output files must all have the same prefixes.")

# from the output files, let us figure out the output prefixes.
if os.sep not in outputs[0]:
    # if no  / separator provided, then outdir is local . directory
    outdir = "."
else:
    # we remove the last file, the remaining prefix is the outdir directory
    outdir = outputs[0].rsplit("/",1)[0]

# ---------------------------------------- manage the input parameters.
params = snakemake.params

qvalue = params.get("qvalue", 0.05)
options = params.get("options", " --keep-dup all ")

# -B --SPMR is to save fragment pileup useful to save into bigwig later on
cmd = "macs3 callpeak -B --SPMR "
cmd += " -t {ip_bam} "
cmd += " -c {control_bam} "
cmd += " -g {params.genome_size} "
cmd += " -n {params.prefix} "
cmd += "--bw {params.bandwidth} "
cmd += " {options} -q {qvalue} "

# is it paired data or not ?
if params.paired:
    cmd += " -f BAMPE "
else:
    cmd += " -f BAM "

# switch between the narrow and broad case
if params.mode == "narrow":
    # First narrow peak calling
    shell(cmd + " --outdir {outdir} {log}")
elif params.mode == "broad":
    # second broad peak calling
    shell(cmd + " --outdir {outdir} --broad --broad-cutoff {params.broad_cutoff} {log}")
else:
    raise ValueError(f"Sequana-wrappers:macs3 error. The 'mode' key in the params section can be only set to 'narrow' or 'broad'; you gave {params.mode}")

