CMD = """\
circular_flag=""
if [ "{params.circular}" = "True" ]; then
    circular_flag="-o"
fi
annotation_flag=""
if [ -n "{params.annotation}" ]; then
    annotation_flag="--annotation-file {params.annotation}"
fi
sequana_coverage \
    --input-file {input.bed} \
    -H {params.high_threshold} \
    -L {params.low_threshold} \
    --clustering-parameter {params.double_threshold} \
    --chunk-size {params.chunksize} \
    --window-gc {params.gc_window_size} \
    --mixture-models {params.mixture_models} \
    --output-directory {params.output_directory} \
    --window-median {params.window_size} \
    --reference-file {input.fasta} \
    $circular_flag \
    $annotation_flag \
    {params.options} \
    >> {log} 2>&1
"""
