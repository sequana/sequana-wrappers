CMD = """\
extra=""
if [ "{params.remove_duplicates}" = "True" ]; then extra="$extra --remove-duplicates"; fi
if [ -n "{params.tmp_directory}" ]; then extra="$extra --tmpdir={params.tmp_directory}"; fi
sambamba markdup {input.bam} {output.bam} {params.options} $extra > {log} 2>&1
"""
