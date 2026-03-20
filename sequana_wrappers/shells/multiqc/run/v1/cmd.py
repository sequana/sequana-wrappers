CMD = """\
module_opts=""
for m in {params.modules}; do
    module_opts="$module_opts -m $m"
done
config_opt=""
if [ -n "{params.config_file}" ]; then
    config_opt="-c {params.config_file}"
fi
multiqc {params.input_directory} \
    --force \
    {params.options} \
    $module_opts \
    $config_opt \
    -o $(dirname {output}) \
    -n $(basename {output}) \
    > {log} 2>&1
"""
