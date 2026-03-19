CMD = """\
options="{params.options}"
if [[ "$options" != *"-F SAF"* ]]; then
    options="$options -t {params.feature} -g {params.attribute}"
fi
featureCounts -T {threads} -s {params.strandness} $options \
    -a {input.gff} -o {output.counts} {input.bam} > {log} 2>&1
"""
