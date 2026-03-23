CMD = """\
alignments="{input.alignments}"
polypolish {params.options} {input.assembly} $alignments 1>{output.fasta} 2>{log}
"""
