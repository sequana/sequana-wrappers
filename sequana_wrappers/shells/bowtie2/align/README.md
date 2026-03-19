# bowtie2/align

Align reads to a reference with Bowtie2, pipe through samtools sort and index.
The index prefix is hardcoded to `reference/genome`.

## Testing

```bash
cd test
make        # run the test (includes bowtie2/build as prerequisite)
make clean  # remove outputs
```

Requires Apptainer/Singularity. The container is pulled automatically on first run.

## Tool versions

| Shell version | Tool | Container |
|---|---|---|
| v1 | bowtie2 2.5.1, samtools | https://zenodo.org/record/8092297/files/bowtie2_2.5.1.img |

## Inputs / outputs

| Name | Description |
|---|---|
| `input.fastq` | One or two FASTQ files (single-end or paired-end) |
| `input.idx` | Bowtie2 index files (triggers build dependency) |
| `output.bam` | Sorted and indexed BAM |

## Params

| Name | Description |
|---|---|
| `params.options` | Extra bowtie2 options |

## Notes

The index prefix is hardcoded to `reference/genome`. Place your reference and
index under `reference/` to match the Sequana mapper pipeline layout.
Single-end and paired-end inputs are detected automatically from the number
of files in `input.fastq`.

## Changelog

- **v1** — initial version: `bowtie2 | samtools sort && samtools index`
