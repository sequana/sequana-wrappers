CMD = """\
db_type=""
case "{output.db}" in
    *.nhr) db_type="nucl" ;;
    *.phr) db_type="prot" ;;
    *) db_type="nucl" ;;
esac
out_name=$(echo "{output.db}" | sed 's/\\.nhr$//' | sed 's/\\.phr$//')
makeblastdb \
    -in {input.fasta} \
    -dbtype $db_type \
    {params.options} \
    -logfile {log} \
    -out $out_name
"""
