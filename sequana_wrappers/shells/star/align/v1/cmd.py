CMD = """\
genomeDir=$(dirname {input.index})
out="{output.bam}"
prefix="${{out%_Aligned.sortedByCoord.out.bam}}"
annot_opt=""
if [ -n "{params.annotation}" ]; then
    annot_opt="--sjdbGTFfile {params.annotation} --sjdbOverhang {params.sjdbOverhang}"
    if [[ "{params.annotation}" == *.gff ]] || [[ "{params.annotation}" == *.gff3 ]]; then
        annot_opt="$annot_opt --sjdbGTFtagExonParentTranscript {params.exonParentTranscript}"
    fi
fi
STAR --genomeDir $genomeDir \
    --twopassMode Basic \
    --twopass1readsN -1 \
    --readFilesIn {input.fastq} \
    --runThreadN {threads} \
    --genomeLoad NoSharedMemory \
    --outSAMtype {params.outSAMtype} \
    --readFilesCommand {params.readFilesCommand} \
    --outFileNamePrefix ${{prefix}}_ \
    $annot_opt \
    {params.options} > {log} 2>&1
"""
