rule rulegraph:
    """Rulegraph 

     Required input:
         - the snakefile filename

     Required output:
         - the output SVG filename 

    Required parameters:
        - mapper: a dictionary mapping each rule to a URL (HTML 
            file or directory). Rules provided in this dictionary will be shown 
            in blue and clickable in the ouptut SVG file.
        - configname: a config file required by the input Snakefile

    """
    input:
        filename = "Snakefile"
    output:
        svg  = "test.svg"
    params:
        mapper = {"rulegraph": "test.html"},
        configname = "config.yaml" 
