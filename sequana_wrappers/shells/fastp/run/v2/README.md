# fastp/run v2

Runs [fastp](https://github.com/OpenGene/fastp) for adapter trimming and quality
filtering. Handles both single-end and paired-end data.

## What changed from v1

`v1` references `{output.r2}` inside the template. Snakemake's formatter resolves
every `{...}` eagerly, so a single-end rule that omits `r2` triggers
`AttributeError: 'OutputFiles' object has no attribute 'r2'` before the bash `if`
block has a chance to run.

`v2` moves the r2 path to `params.out2`:

- single-end rules set `params.out2 = ""` and never touch `output.r2`
- paired-end rules set `params.out2 = lambda wc, output: output.r2` (or the path
  directly) and declare `output.r2` as usual

The template never references `{output.r2}`, so single-end rules no longer break.

## Rule interface

### input

| Name    | Description                                    |
|---------|------------------------------------------------|
| `fastq` | List of input FastQ files (1 for SE, 2 for PE) |

### output

| Name   | Required    | Description                      |
|--------|-------------|----------------------------------|
| `r1`   | always      | Trimmed R1 (or single-end) FastQ |
| `r2`   | paired only | Trimmed R2 FastQ                 |
| `html` | always      | fastp HTML report                |
| `json` | always      | fastp JSON report (for MultiQC)  |

### params

| Name       | Description                                                    |
|------------|----------------------------------------------------------------|
| `options`  | Extra fastp CLI options (string)                               |
| `adapters` | Adapter options string (may be `""`)                           |
| `out2`     | R2 output path (PE) or `""` (SE). Lambda form recommended.     |

### log

A single log file. stdout and stderr are both redirected (`> {log} 2>&1`).

## Example rule (single-end)

```python
rule fastp:
    input:
        fastq=manager.getrawdata()
    output:
        r1="{sample}/fastp/{sample}.fastq.gz",
        html="{sample}/fastp/fastp_{sample}.html",
        json="{sample}/fastp/fastp_{sample}.json",
    log:
        "logs/fastp/{sample}.log"
    params:
        options=config["fastp"]["options"],
        adapters=config["fastp"]["adapters"],
        out2=""
    threads:
        config["fastp"].get("threads", 4)
    container:
        config["apptainers"]["fastp"]
    shell:
        manager.get_shell("fastp/run", "v2")
```

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
        out2=lambda wildcards, output: output.r2
    threads:
        config["fastp"].get("threads", 4)
    container:
        config["apptainers"]["fastp"]
    shell:
        manager.get_shell("fastp/run", "v2")
```
