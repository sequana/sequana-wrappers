# fastp/run v1

Runs [fastp](https://github.com/OpenGene/fastp) for adapter trimming and quality filtering.
Handles both single-end and paired-end data by counting the number of input files.

## Rule interface

The Snakemake rule using this shell **must** declare the following named fields:

### input

| Name     | Description                                    |
|----------|------------------------------------------------|
| `fastq`  | List of input FastQ files (1 for SE, 2 for PE) |

### output

| Name     | Required      | Description                        |
|----------|---------------|------------------------------------|
| `r1`     | always        | Trimmed R1 (or single-end) FastQ   |
| `r2`     | paired only   | Trimmed R2 FastQ                   |
| `html`   | always        | fastp HTML report                  |
| `json`   | always        | fastp JSON report (for MultiQC)    |

### params

| Name       | Description                          |
|------------|--------------------------------------|
| `options`  | Extra fastp CLI options (string)     |
| `adapters` | Adapter options string (may be `""`) |

### log

A single log file. Both stdout and stderr are redirected to it (`> {log} 2>&1`).

## Example rule (paired-end)

```python
rule fastp:
    input:
        fastq=manager.getrawdata()
    output:
        r1="{sample}/fastp/{sample}_R1_.fastp.fastq.gz",
        r2="{sample}/fastp/{sample}_R2_.fastp.fastq.gz",
        html="{sample}/fastp/fastp_{sample}.html",
        json="{sample}/fastp/fastp_{sample}.json",
    log:
        "{sample}/fastp/{sample}.log"
    params:
        options=config["fastp"]["options"],
        adapters=config["fastp"]["adapters"],
    threads:
        config["fastp"].get("threads", 4)
    container:
        config["apptainers"]["fastp"]
    shell:
        manager.get_shell("fastp/run", "v1")
```

## Notes

- `r2` must **not** be declared for single-end data — the bash `if` block checks the number
  of input files at runtime, so only `r1` is needed in the single-end rule definition.
- `params.adapters` is inserted verbatim between `{params.options}` and the read arguments;
  pass `""` to disable explicit adapter specification (fastp will auto-detect adapters).
