CMD = """\
alignments="{input.alignments}"
polypolish polish {params.options} {input.assembly} $alignments 1>{output.fasta} 2>{log}
"""
