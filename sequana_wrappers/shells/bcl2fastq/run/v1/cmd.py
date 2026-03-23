CMD = """\
bcl2fastq -p {threads} \
    --barcode-mismatches {params.barcode_mismatch} \
    --runfolder-dir {params.indir} \
    --intensities-dir {params.indir}/Data/Intensities \
    --sample-sheet {input.samplesheet} \
    --output-dir {params.outdir} \
    {params.options} > {log} 2>&1
"""
