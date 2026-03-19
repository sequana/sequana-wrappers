# bowtie2/build

Build a Bowtie2 index from a reference FASTA. The index is always written to
`reference/genome` (prefix hardcoded to match the mapper pipeline convention).

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
| v1 | bowtie2 2.5.1 | https://zenodo.org/record/8092297/files/bowtie2_2.5.1.img |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.reference` | Reference FASTA (expected at `reference/genome.fasta`) |
| `output` | Six Bowtie2 index files (`reference/genome.{1..4}.bt2`, `.rev.{1,2}.bt2`) |

## Params

| Name | Description |
|---|---|
| `params.options` | Extra bowtie2-build options |

## Notes

The output prefix is hardcoded to `reference/genome` to match the Sequana
mapper pipeline layout. The reference FASTA must therefore be placed at
`reference/genome.fasta`.

## Changelog

- **v1** — initial version: `bowtie2-build`
