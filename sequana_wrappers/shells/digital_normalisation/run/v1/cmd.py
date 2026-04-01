CMD = """\
fastq_files=({input.fastq})
if [ ${{#fastq_files[@]}} -eq 2 ]; then
    r1="${{fastq_files[0]%.gz}}"
    r2="${{fastq_files[1]%.gz}}"
    [[ "${{fastq_files[0]}}" == *.gz ]] && unpigz -p {threads} -fk "${{fastq_files[0]}}"
    [[ "${{fastq_files[1]}}" == *.gz ]] && unpigz -p {threads} -fk "${{fastq_files[1]}}"
    prefix=$(basename "$r1" .fastq)
    interleave-reads.py "$r1" "$r2" --output ${{prefix}}.pe > {log} 2>&1
    normalize-by-median.py --paired --ksize {params.ksize} \
        --cutoff {params.cutoff} -M {params.m} {params.options} \
        --savegraph {params.tmp_graph_name} ${{prefix}}.pe \
        --output ${{prefix}}.pe.keep >> {log} 2>&1
    filter-abund.py --threads {threads} -V {params.tmp_graph_name} \
        ${{prefix}}.pe.keep --output ${{prefix}}.pe.filter >> {log} 2>&1
    extract-paired-reads.py ${{prefix}}.pe.filter \
        --output-paired ${{prefix}}.pe \
        --output-single ${{prefix}}.orphans >> {log} 2>&1
    split-paired-reads.py ${{prefix}}.pe \
        -1 {output.r1} -2 {output.r2} >> {log} 2>&1
    rm -f ${{prefix}}.pe ${{prefix}}.pe.filter ${{prefix}}.pe.keep ${{prefix}}.orphans {params.tmp_graph_name}
    [[ "${{fastq_files[0]}}" == *.gz ]] && rm -f "$r1" "$r2"
else
    f="${{fastq_files[0]}}"
    [[ "$f" == *.gz ]] && unpigz -p {threads} -fk "$f" && f="${{f%.gz}}"
    prefix=$(basename "$f" .fastq)
    normalize-by-median.py --ksize {params.ksize} \
        --cutoff {params.cutoff} -M {params.m} {params.options} \
        --savegraph {params.tmp_graph_name} "$f" \
        --output ${{prefix}}.se.keep >> {log} 2>&1
    filter-abund.py --threads {threads} -V {params.tmp_graph_name} \
        ${{prefix}}.se.keep --output ${{prefix}}.se.filter >> {log} 2>&1
    mv ${{prefix}}.se.filter {output.r1}
    rm -f ${{prefix}}.se.keep {params.tmp_graph_name}
    [[ "${{fastq_files[0]}}" == *.gz ]] && rm -f "$f"
fi
"""
