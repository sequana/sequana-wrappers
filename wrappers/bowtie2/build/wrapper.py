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
from pathlib import Path
from snakemake.shell import shell

options = snakemake.params.get("options", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
done = snakemake.output.done
reference = snakemake.input.reference

#extract reference parent path and genome name to set indexbase
genome_directory = Path(reference).parent
genome_name = Path(reference).name.rsplit(".")[0]


# now let us figure out the version of bowtie2.
try:
    p = subprocess.Popen(['bowtie2'], stderr=subprocess.PIPE)
    p.wait()
    stderr = p.stderr.read().decode().split("\n")
    hits = [line for line in stderr if "version" in line and 'Bowtie' in line]
    bowtie2_version = "_" + hits[0].split("version")[1].split()[0].strip()
except Exception:  # various type of exception may occur here
     logger.warning(f"Could not determine bowtie2 version. Index will be stored in {genome_directory}/bowtie2/    ")
     bowtie2_version = ""

index_base = f"{genome_directory}/bowtie2{bowtie2_version}/{genome_name}"




shell("bowtie2-build --threads {snakemake.threads} {options} {snakemake.input.reference} {indexbase} {log}")
shell("samtools faidx {snakemake.input.reference}")
# this is required to make the rule complete
shell("touch {done})

