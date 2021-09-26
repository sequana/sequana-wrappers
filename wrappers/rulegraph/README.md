
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

# Example

    rulegraph_params_mapper = {"a_given_rule": "its_html.html"}

    rule rulegraph:
        input: manager.snakefile
        output:
            svg = ".sequana/rulegraph.svg"
        params:
            mapper = rulegraph_params_mapper,
            configname = "config.yaml"
        wrapper:
            "main/wrappers/rulegraph"

# Requirements

- dot
