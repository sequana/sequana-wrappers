# The Sequana Wrapper Repository

[![Tests wrappers](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml/badge.svg)](https://github.com/sequana/sequana-wrappers/actions/workflows/main.yml)
[![Tests shells](https://github.com/sequana/sequana-wrappers/actions/workflows/shells.yml/badge.svg)](https://github.com/sequana/sequana-wrappers/actions/workflows/shells.yml)
[![Tests](http://joss.theoj.org/papers/10.21105/joss.00352/status.svg)](http://joss.theoj.org/papers/10.21105/joss.00352)

|||
| --- | --- |
| Overview | Shell command library and Snakemake wrappers for Sequana pipelines |
| Status | Production (`wrappers/` — maintenance only) / Active (`shells/`) |
| Issues | Please fill a report on [github/sequana/sequana-wrappers](https://github.com/sequana/sequana/issues) |
| Python version | Python 3.8+ |
| Citation | Cokelaer et al, (2017), ‘Sequana’: a Set of Snakemake NGS pipelines, Journal of Open Source Software, 2(16), 352, [doi:10.21105/joss.00352](http://www.doi2bib.org/bib/10.21105%2Fjoss.00352) |

## Status and roadmap

This repository contains two independent mechanisms for providing bioinformatics
tool commands to Sequana pipelines:

- **`wrappers/`** — the original Snakemake wrapper system (Python scripts +
  conda `environment.yaml`).  **This tree is now in maintenance mode.**
  No new wrappers will be added.  Bug fixes will still be accepted, but all
  new development happens in `shells/`.  See [the rationale below](#the-shells-directory--rationale-and-design)
  for the full explanation of why.

- **`sequana_wrappers/shells/`** — the new shell command library.  Versioned
  shell strings that work with `container:` + `shell:` rules, with no Python
  inside the container.  This is the active development track and the
  recommended approach for all new Sequana pipelines.

## Quick start — shells (recommended)

Install the package:

```bash
pip install sequana_wrappers
```

Use in a Snakemake pipeline via `sequana_pipetools`:

```python
# In your pipeline’s .rules file — manager is a PipelineManager instance
rule minimap2:
    input:   ...
    output:  "{sample}/{sample}.sorted.bam"
    container: "https://zenodo.org/record/7987999/files/samtools_1.17_minimap2_2.24.0.img"
    shell:   manager.get_shell("minimap2/align", "v1")
```

## Quick start — wrappers (legacy)

```bash
snakemake --wrapper-prefix https://github.com/sequana/sequana-wrappers
```

or with a local copy:

```bash
git clone git@github.com:sequana/sequana-wrappers.git sequana_wrappers
snakemake --wrapper-prefix git+file:///home/user/sequana_wrappers
```

If the environment variable `SEQUANA_WRAPPERS` is set to
`git+file:///home/user/sequana_wrappers`, all pipelines will automatically use
it as the `--wrapper-prefix`.

# The `shells/` directory — rationale and design

## Background

Sequana pipelines use two mechanisms to provide bioinformatics tools:

- **Wrappers** (`wrappers/`) — Python scripts (`wrapper.py`) plus a conda
  `environment.yaml`.  Snakemake fetches and executes them via the `wrapper:`
  rule directive.
- **Containers** — Apptainer/Singularity images (hosted on Zenodo/Damona)
  referenced by the `container:` rule directive.

## The problem: `wrapper:` + `container:` are incompatible

When a Snakemake rule combines both `wrapper:` and `container:`, and the
pipeline is run with `--use-singularity` / `--apptainer-prefix`, Snakemake v7
executes the wrapper Python script **inside the container** using the
container's own Python binary.

Snakemake does not require Snakemake to be installed inside the container — it
bind-mounts its own `site-packages` from the host at `/mnt/snakemake`.
However, the **container's Python binary** must be ABI-compatible with the
host Python.  Old bioconda/Damona images (Python 3.8) fail when the host runs
Python 3.10 because C-extension `.so` files compiled for 3.10 cannot be loaded
by Python 3.8:

```
ImportError: /mnt/snakemake/...so: cannot open shared object
```

A further sign of inverted concerns: **Damona had to ship Python inside
tool containers** (bwa, samtools, …) specifically to satisfy the wrapper
mechanism.  A container for `bwa` should contain `bwa`, not a Python runtime
serving the pipeline framework.

## Options considered

### Option 1 — Update container images to Python 3.10
Rebuild all Damona images with a Python version matching the host.  Wrappers
would then work with containers as originally intended.

*Rejected as the primary fix*: containers would still carry Python only to
serve the framework; images must be rebuilt every time the host Python is
upgraded; the inverted concern is not resolved.

### Option 2 — Remove `container:` from wrapper rules (short-term workaround)
Wrapper rules run on the host (or via `--use-conda`); only pure `shell:` rules
keep their `container:` directive.

*Used as a temporary workaround* in sequana_mapper while the shell library
was being designed.  Downside: Apptainer only covers a subset of rules.

### Option 3 — Drop wrappers, inline `shell:` in each pipeline (simplest)
Replace every wrapper with a hand-written `shell:` block inside the pipeline
rule.  No shared library.

*Rejected*: duplicates logic across all pipelines; maintenance burden; loses
the reusability benefit of this repository entirely.

### Option 4 — Shell command library in `shells/` ✓ (chosen)
Return to the spirit of the early sequana approach: define reusable, versioned
shell command strings here, alongside the existing `wrappers/`.  Pipelines
import these strings and use them in `shell:` + `container:` rules.

**Why this wins:**

| Property | Wrappers | Shell library |
|---|---|---|
| Reusable logic | Yes (Python) | Yes (string) |
| Python in container | Required | **Not needed** |
| Git tag checkout at run time | Yes | **No** |
| Damona images lean | No | **Yes** |
| Works with `--use-conda` | Yes | No |
| Apptainer compatible | Only if Python ABI matches | **Always** |
| Backward compatible | — | Yes (`wrappers/` kept) |

The `wrappers/` tree is kept untouched for full backward compatibility.

## Design

### Repository layout

```
sequana-wrappers/
├── wrappers/                        # existing — kept for backward compat
│   ├── bwa/align/wrapper.py
│   └── ...
└── sequana_wrappers/shells/         # Python package — container-first shell library
    ├── __init__.py
    ├── bwa/
    │   ├── __init__.py
    │   ├── align/
    │   │   ├── __init__.py
    │   │   └── v1/
    │   │       ├── __init__.py
    │   │       └── cmd.py           # frozen at release v1
    │   └── build/
    │       ├── __init__.py
    │       └── v1/cmd.py
    ├── bamtools/
    │   └── stats/
    │       ├── __init__.py
    │       └── v1/cmd.py
    └── ...
```

### Versioning convention

Every shell script is named `cmd.py` and lives inside a named version
subdirectory.  The structure is:

```
sequana_wrappers/shells/<tool>/<command>/<version>/cmd.py
```

Valid version names are:

- **`vN`** (e.g. `v1`, `v2`) — frozen, reproducible snapshots.  Once created,
  these files are **never edited**.
- **`dev`** — work-in-progress version used during active development.
  A `dev/` directory is created when new work begins on a command and removed
  (or renamed to `vN`) at release time.  **No `dev/` directories exist in
  released versions of this package.**

Every shell script is named `cmd.py`; the tool and command are encoded
entirely in the directory path.  This makes future deeper nesting
(e.g. `shells/bamtools/stats/paired/v1/cmd.py`) natural without any
changes to the `get_shell` API.

There is **no silent fallback** between versions: requesting a version that
does not exist raises an explicit error.

### Version axes

Two version axes are completely independent:

| Axis | What it pins | Where it lives |
|---|---|---|
| Tool binary | e.g. `bwa 0.7.17` | Container image (Damona / Zenodo) |
| Shell command | e.g. `v1` | hardcoded per rule in the pipeline |

Each pipeline rule hardcodes its own shell command version independently —
one rule can use `v1` while another uses `v2` if only that command changed.

### Shell file format

Each `cmd.py` exports a single `CMD` string using Snakemake's standard
`{input}`, `{output}`, `{params}`, `{threads}`, `{log}`, and `{wildcards}`
placeholders:

```python
# sequana_wrappers/shells/bwa/align/v1/cmd.py
CMD = """\
mkdir -p {params.tmp_directory}
(bwa mem -t {threads} {params.options} {input.reference} {input.fastq} \
 | sambamba view -t {threads} -S -f bam -o /dev/stdout /dev/stdin \
 | sambamba sort /dev/stdin -o {output.sorted} -t {threads} \
   --tmpdir={params.tmp_directory}) \
> {log} 2>&1
"""
```

### Usage in a pipeline rule

`get_shell` is available as a method on the `PipelineManager` instance
(already present in every pipeline rules file) — no extra import needed.
The version is hardcoded per rule:

```python
rule bwa:
    input:   ...
    output:  sorted="{sample}/{sample}.sorted.bam"
    log:     "{sample}/bwa/{sample}.log"
    params:  options=config["bwa"]["options"],
             tmp_directory=config["bwa"]["tmp_directory"]
    threads: 2
    container: config['apptainers']['bwa']
    shell:   manager.get_shell("bwa/align", "v1")
```

Use `"dev"` during development before a versioned snapshot exists:

```python
    shell:   manager.get_shell("bwa/align", "dev")
```

The container contains **only the tool binaries** — no Python, no Snakemake.

### Adding or updating a shell command

**During development:**

1. Create `sequana_wrappers/shells/<tool>/<command>/dev/` with `__init__.py`
   and `cmd.py`.
2. Use `manager.get_shell("<tool>/<command>", "dev")` in the pipeline rule.
3. Test against the relevant Damona container.

**At release time:**

```bash
VERSION=v2
mkdir -p sequana_wrappers/shells/<tool>/<command>/${VERSION}
touch sequana_wrappers/shells/<tool>/<command>/${VERSION}/__init__.py
cp sequana_wrappers/shells/<tool>/<command>/dev/cmd.py \
   sequana_wrappers/shells/<tool>/<command>/${VERSION}/cmd.py
rm -rf sequana_wrappers/shells/<tool>/<command>/dev
```

Then update the pipeline rule to `manager.get_shell("<tool>/<command>", "v2")`
and bump `version` in `pyproject.toml`.  No git tag required — the directory
**is** the version.

---

# Notes for developers

## Overview

> **`wrappers/` is in maintenance mode.**  Bug fixes are welcome; new wrappers
> are not accepted.  All new tool commands should be added to
> `sequana_wrappers/shells/` instead.  See the
> [shells rationale](#the-shells-directory--rationale-and-design) for context.

The `wrappers/` directory contains the legacy wrappers. Each sub-directory is dedicated to
a wrapper that is related to a given software/application. A sub directory may have several wrappers (e.g., bwa has a sub directory related to the indexing, and a sub directory related to mapping).

Here is an example of a wrapper tree structure:

    fastqc
    ├── environment.yaml
    ├── README.md
    ├── test
    │   ├── README.md
    │   ├── Snakefile
    │   ├── test_R1_.fastq
    │   └── test_R2_.fastq
    └── wrapper.py

Note that some software may have several sub wrappers (see the bowtie1 wrapper for instance).

A wrapper directory must contain a file called **wrapper.py** where the
developers must provide the core of the wrapper. There is no specific
instructions here except to write good code as much as possible (with comments).

A wrapper directory should have a **test** directory for continuous integration
with a **Snakefile** to be tested and possibly data file **Do not add large files here**. A
**README.md** should be added to explain the origin of the test data files.
Finally, include your tests in the main [**test.py**](test.py) file
of the root of the repository (not the wrapper itself).

For testing purposes, you should also add a file called **environment.yaml**
to tell what are the required packages to be installed for the test (and wrapper)
to work.

Finally, for the documentation, we ask the developer to create a **README.md** file
described here below.

To test your new wrapper (called *example* here), type:

   pytest test.py -k test_example

## The config file

If required in a wrapper, parameters must be defined in a **config.yaml** file.
Similarly for threading. Consider the following pointswhen writting a wrapper:

- The thread paramter should also be a parameter in config file.
- the **params** section should contain a key called **options** also define in the config file.
- keys or parameters related to directories and files should use the *_directory* or *_file* suffices. This is for
  Sequanix application to automatically recognised those options with a dedicated widget.

Consider this example:

    rule falco:
        input: whatever
        output: whatever
        log:
            "samples/{sample}/falco.log"
        threads:
            config['falco']['threads']
        params:
            options=config['falco']['options'],
            wkdir=config['falco']['working_directory']
        wrapper:
            "falco/wrappers/falco"

You config file will look like:

    falco:
        threads: 4
        options="--verbose"
        working_directory: 'here'


## Naming arguments of the different sections

In all sections (e.g., input, params), if there is only one input, no need to name it, otherwise, please do.

    rule example1:
        input:
            "test.bam"
        output:
            "test.sorted.bam"
        ...

but:

    rule example1:
        input:
            "test.bam"
        output:
            bam="test.sorted.bam"
            bai="test.sorted.bam.bai"
        ...




## Documentation

Each wrapper should have a dedicated documentation explaining the input/output with a usage example. It should also document the expected configuration file.  The file must be formatted in markdown. It must contain a **Documentation** and **Example** sub sections. If a **Configuration** section is found, it is also added to the documentation. This **README.md**  file will be rendered automatically via a Sequana sphinx plugin. Consider the [fastqc](wrappers/fastqc/README.md) directory for a workable example rendered [here](https://sequana.readthedocs.io/en/master/wrappers.html#fastqc).



## Faqs

### adding a new wrapper in practice

In ./wrappers, add a new wrapper. Copy the existing fastqc wrapper for instance.
Edit the wrapper.py and design a test/Snakefile example for testing. Since you are
a developer, you are problaby developping in a dedicated branch. Let us call it **dev**.

In the test/Snakefile, you should switch from the **main** to the **dev** in the wrapper path:

    wrapper:
        "dev/wrappers/my_new_wrapper"

In order to test your Snakefile, you first need to commit the wrapper.py. Then, execute the Snakefile:

    snakemake -s Snakefile  -j 1 --wrapper-prefix git+file:///YOURPATH/sequana-wrappers/ -f -p

If it fails, edit and commit your wrapper.py and execute again until your Snakefile and wrappers are functional.

Once done, switch back the wrapper path to the **main** branch:

    wrapper:
        "main/wrappers/my_new_wrapper"

Time to include the new wrapper in the continous integration. Go to the root of sequana-wrappers and add a functional test to the end of test.py. Then, test it:

    pytest test.py -k my_new_wrapper -v

You are ready to push and create a pull-requests



## Tagging

You may consider to add a tag. Our convention is to use a tag with the YEAR.MONTH.DAY
where day and month do not include extra zeros. So you would have e.g.:

    v23.11.11
    v23.2.2

but not v23.02.02

We use annotated tags so the command is e.g.:

    git tag -a v23.11.11

and then push it to github:

    git push origin main v23.11.11


