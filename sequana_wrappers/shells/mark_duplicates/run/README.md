# mark_duplicates/run

Mark duplicate reads in a coordinate-sorted BAM file using Picard MarkDuplicates,
then index the output with samtools.

## Usage

```python
from sequana_wrappers import get_shell

rule mark_duplicates:
    input:
        bam="sample/sample.bam"
    output:
        bam="sample/sample.markdup.bam",
        metrics="sample/sample.markdup.metrics"
    log:
        "logs/mark_duplicates.log"
    params:
        options="",
        remove_dup="false",
        tmpdir="tmp"
    threads: 4
    container:
        "https://zenodo.org/record/18257162/files/sequana_tools_26.1.14.img"
    shell:
        get_shell("mark_duplicates/run", "v1")
```

## Parameters

| Parameter | Description |
|---|---|
| `params.remove_dup` | Remove duplicates instead of flagging them (`"true"` or `"false"`) |
| `params.tmpdir` | Temporary directory for Picard |
| `params.options` | Extra Picard MarkDuplicates options |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.bam` | Coordinate-sorted BAM file with read group (`@RG`) information |
| `output.bam` | BAM file with duplicates flagged (or removed) |
| `output.metrics` | Picard duplication metrics file |

The output BAM is indexed (`.bai`) automatically.

## Testing

```bash
cd test
make        # run the test
make clean  # remove outputs
```

Requires Apptainer/Singularity. The container is pulled automatically on first run.

The input BAM (`test/sample/sample.bam`) must contain read group information (`@RG` header).

## Tool versions

| Shell version | Tool | Container |
|---|---|---|
| v1 | Picard 3.4.0, samtools | https://zenodo.org/record/18257162/files/sequana_tools_26.1.14.img |

## Changelog

- **v1** — initial version: `picard MarkDuplicates` + `samtools index`
