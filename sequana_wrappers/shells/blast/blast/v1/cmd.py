CMD = """\
{params.blast_type} \
    -query {input.query} \
    -db {params.db} \
    -num_threads {threads} \
    -outfmt '{params.outfmt}' \
    -evalue {params.evalue} \
    {params.options} \
    -out {output} 2> {log}
"""
