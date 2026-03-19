# bamtools/stats

Compute BAM statistics (total reads, mapped, duplicates, insert size, etc.)
using bamtools stats.

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
| v1 | bamtools 2.5.2 | https://zenodo.org/record/10211475/files/bamtools_2.5.2.img |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.bam` | Input BAM file |
| `output` | Text file with statistics |

## Changelog

- **v1** — initial version: `bamtools stats -in {input.bam} -insert`
