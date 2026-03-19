# bwa/build

Index a reference FASTA with BWA and create a samtools `.fai` index.

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
| v1 | bwa 0.7.17, samtools | https://zenodo.org/record/7341710/files/sequana_tools_0.14.5.img |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.reference` | Reference FASTA |
| `output.bwa_bwt` | BWA index (`.bwt`) |
| `output.fai` | Samtools FASTA index (`.fai`) |

## Params

| Name | Description |
|---|---|
| `params.index_algorithm` | `is` (short genomes) or `bwtsw` (long genomes) |
| `params.options` | Extra bwa index options |

## Changelog

- **v1** — initial version: `bwa index` + `samtools faidx`
