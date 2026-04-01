CMD = """\
wkdir=$(dirname {output.done})
gff_opt=""
if [ -n "{input.gff}" ]; then
    gffread {input.gff} -g {input.fasta} -w $wkdir/salmon_transcript.fa > {log} 2>&1
    fasta_in=$wkdir/salmon_transcript.fa
else
    fasta_in={input.fasta}
fi
salmon index -p {threads} -t $fasta_in -i $wkdir {params.options} >> {log} 2>&1
touch {output.done}
"""
