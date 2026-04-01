CMD = """\
out="{output[0]}"
out="${{out%/}}"
out_path=$(dirname "$out")
out_name=$(basename "$out")
path_opt=""
if [ "$out_path" != "." ]; then
    path_opt="--out_path $out_path"
fi
download_opt=""
if [ -n "{params.downloads_path}" ]; then
    download_opt="--download_path {params.downloads_path}"
fi
busco --in {input} --out $out_name --force \
    $path_opt \
    --cpu {threads} \
    --mode {params.mode} \
    --lineage {params.lineage} \
    $download_opt \
    {params.options} > {log} 2>&1
"""
