rule rulegraph:
    """Creates a rulegraph showing your pipeline dependencies

     Required input:
         - the snakefile filename

     Required output:
         - **svg**: the output SVG filename

     Required parameters:
        - **mapper**: a dictionary mapping each rule to a URL (HTML
            file or directory). Rules provided in this dictionary will be
            shown in blue and clickable in the ouptut SVG file.
        - **configname**: a config file required by the input Snakefile

    """
    input: manager.snakefile
    output:
        svg = ".sequana/rulegraph.svg"
    params:
        mapper = __rulegraph__mapper,
        configname = "config.yaml"
    wrapper:
        "main/wrappers/rulegraph"

