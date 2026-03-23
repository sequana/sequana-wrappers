CMD = """\
dsrc d -s -t{threads} {params.options} {input} | pigz -p {threads} > {output}
pigz -p {threads} --test {output}
rm -f {input}
"""
