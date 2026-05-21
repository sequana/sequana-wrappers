CMD = """\
fastq_files=({input.fastq})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    input_opt="-1 ${{fastq_files[0]}} -2 ${{fastq_files[1]}}"
else
    input_opt="-s ${{fastq_files[0]}}"
fi
preset_opt=""
if [ -n "{params.preset}" ]; then preset_opt="--{params.preset}"; fi
outdir=$(dirname {output.contigs})
spades.py ${{preset_opt}} \
    -k {params.k} \
    --memory {params.memory} \
    --threads {threads} \
    ${{input_opt}} \
    -o ${{outdir}} \
    {params.options} > {log} 2>&1
cp ${{outdir}}/contigs.fasta {output.contigs}
cp ${{outdir}}/scaffolds.fasta {output.scaffolds}
"""
