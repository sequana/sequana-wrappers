import os
import tempfile
from pathlib import Path

import pytest

snakemake = pytest.importorskip("snakemake")

from sequana_wrappers.snippets.rulegraph.run.v1.code import execute


SIMPLE_SNAKEFILE = """\
rule all:
    input: "output.txt"

rule make_output:
    output: "output.txt"
    shell: "touch {output}"
"""


def test_execute_creates_dot_file():
    origdir = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdir:
        snakefile = Path(tmpdir) / "Snakefile"
        snakefile.write_text(SIMPLE_SNAKEFILE)
        output_dot = Path(tmpdir) / "rulegraph.dot"

        os.chdir(tmpdir)
        try:
            execute(
                input=[str(snakefile)],
                output=[str(output_dot)],
                params={},
            )
        finally:
            os.chdir(origdir)

        assert output_dot.exists(), f"Expected output dot file at {output_dot}"
