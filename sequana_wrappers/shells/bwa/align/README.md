# bwa/align

Align reads to a reference with BWA-MEM, then sort and convert to BAM with sambamba.

## Testing

```bash
cd test
make        # run the test (includes bwa/build as prerequisite)
make clean  # remove outputs
```

Requires Apptainer/Singularity. The container is pulled automatically on first run.

## Tool versions

| Shell version | Tool | Container |
|---|---|---|
| v1 | bwa 0.7.17, sambamba | https://zenodo.org/record/7341710/files/sequana_tools_0.14.5.img |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.fastq` | One or two FASTQ files (single-end or paired-end) |
| `input.reference` | Reference FASTA |
| `input.bwa_bwt` | BWA index `.bwt` file (triggers indexing dependency) |
| `input.fai` | Samtools `.fai` index (triggers faidx dependency) |
| `output.sorted` | Sorted BAM |

## Params

| Name | Description |
|---|---|
| `params.options` | Extra bwa mem options (e.g. `-T 30 -M`) |
| `params.tmp_directory` | Temporary directory for sambamba sort |

## Changelog

- **v1** — initial version: `bwa mem | sambamba view | sambamba sort`
