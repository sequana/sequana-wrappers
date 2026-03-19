CMD = """\
module_opts=""
for m in {params.modules}; do
    module_opts="$module_opts -m $m"
done
multiqc {params.input_directory} \
    --force \
    {params.options} \
    $module_opts \
    -o $(dirname {output}) \
    -n $(basename {output}) \
    > {log} 2>&1
"""
