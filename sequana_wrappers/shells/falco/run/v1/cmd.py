CMD = """\
falco -t {threads} --outdir {params.working_directory} {input.fastq} {params.options} > {log} 2>&1
"""
