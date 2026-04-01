CMD = """\
pbunzip2 -p {threads} --test {input[0]} && \
pbunzip2 -p {threads} {input[0]} | pigz -p {threads} > {output[0]} && \
pigz -p {threads} --test {output[0]} && \
rm -f {input[0]}
"""
