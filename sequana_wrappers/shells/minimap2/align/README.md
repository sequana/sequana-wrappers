# minimap2/align

Map reads to a reference with minimap2, pipe through samtools sort.

## Testing

```bash
cd test
make        # run the test
make clean  # remove outputs
```

Requires Apptainer/Singularity. The container is pulled automatically on first run.

## Tool versions

| Shell version | Tool | Container |
|---|---|---|
| v1 | minimap2 2.24.0, samtools | https://zenodo.org/record/7817800/files/minimap2_2.24.0.img |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.fastq` | One or two FASTQ files (single-end or paired-end) |
| `input.reference` | Reference FASTA |
| `output` | Sorted BAM |

## Params

| Name | Description |
|---|---|
| `params.options` | Extra minimap2 options (e.g. `-x map-ont`) |

## Changelog

- **v1** — initial version: `minimap2 | samtools sort`
