# Documentation

Creates a rulegraph showing your pipeline dependencies

Input section:

- the snakefile filename

Output section:

- **svg**: the output SVG filename

Params section:

- **mapper**: a dictionary mapping each rule to a URL (HTML
  file or directory). Rules provided in this dictionary will be
  shown in blue and clickable in the ouptut SVG file.
- **configname**: a config file required by the input Snakefile
- **required_local_files**: a list of required files and directories next to
  the Snakefile that are required to run the pipeline.

# Example

    rulegraph_params_mapper = {"a_given_rule": "its_html.html"}

    rule rulegraph:
        input: 
            manager.snakefile
        output:
            svg = ".sequana/rulegraph.svg"
        params:
            configname = "config.yaml",
            mapper = rulegraph_params_mapper,
            required_local_files = ['rules/',]
        wrapper:
            "main/wrappers/rulegraph"

# Requirements

- dot
