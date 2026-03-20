import os
from pathlib import Path

import pytest

from sequana_wrappers import get_run


SNAKEFILE_CONTENT = """\
rule all:
    input: "output.txt"

rule create:
    output: "output.txt"
    shell: "touch {output}"
"""


def test_get_run_rulegraph_returns_callable():
    fn = get_run("rulegraph/run", "v1")
    assert callable(fn)


def test_rulegraph_execute(tmp_path):
    snakefile = tmp_path / "Snakefile"
    snakefile.write_text(SNAKEFILE_CONTENT)
    output_dot = tmp_path / "rulegraph.dot"

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        execute = get_run("rulegraph/run", "v1")
        execute(
            input=[str(snakefile)],
            output=[str(output_dot)],
            params={},
        )
    finally:
        os.chdir(cwd)

    assert output_dot.exists(), "rulegraph dot file was not produced"
