CMD = """\
mkdir -p {params.wkdir}
genome_length=$(awk '/^[^>]/ {{s+=length($0)}} END {{print s}}' {input.fasta})
nbases=$(awk -v gl="$genome_length" 'BEGIN {{n=int(log(gl)/log(2)/2-1); n=(n<14?n:14); print n}}')
annot_opt=""
if [ -n "{params.annotation}" ]; then
    annot_opt="--sjdbGTFfile {params.annotation} --sjdbOverhang {params.sjdbOverhang}"
    if [[ "{params.annotation}" == *.gff ]] || [[ "{params.annotation}" == *.gff3 ]]; then
        annot_opt="$annot_opt --sjdbGTFtagExonParentTranscript {params.exonParentTranscript}"
    fi
fi
STAR --runMode genomeGenerate \
    --genomeFastaFiles {input.fasta} \
    --genomeDir {params.wkdir} \
    --runThreadN {threads} \
    --genomeSAindexNbases $nbases \
    $annot_opt \
    {params.options} > {log} 2>&1
touch {output}
"""
