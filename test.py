# Testing code copied from snakemake-wrappers
# https://github.com/snakemake/snakemake-wrappers/blob/master/test.py

import subprocess
import os
import tempfile
import shutil

import pytest


DIFF_MASTER = os.environ.get("DIFF_MASTER", "false") == "true"
DIFF_LAST_COMMIT = os.environ.get("DIFF_LAST_COMMIT", "false") == "true"

if DIFF_MASTER or DIFF_LAST_COMMIT:
    compare = "HEAD^" if DIFF_LAST_COMMIT else "origin/main"

    # check if wrapper is modified compared to master
    DIFF_FILES = subprocess.check_output(["git", "diff", compare, "--name-only"]).decode().split("\n")


# what to skip on github or locally
WRAPPERS_TO_SKIP_ON_GA = {"wrappers/bcl2fastq"}


class Skipped(Exception):
    pass


skip_if_not_modified = pytest.mark.xfail(raises=Skipped)

command = ["snakemake", "--cores", "1", "-F", "--use-conda"]


# a copy wrapper function
def copy_wrapper(wrapper, dst):
    copy = lambda pth, src: shutil.copy(os.path.join(pth, src), os.path.join(dst, pth))
    success = False
    for ext in ("py", "R", "Rmd"):
        script = "wrapper." + ext
        if os.path.exists(os.path.join(wrapper, script)):
            os.makedirs(os.path.join(dst, wrapper), exist_ok=True)
            copy(wrapper, script)
            success = True
            break
    assert success, "No wrapper script found for {}".format(wrapper)
    if os.path.exists(wrapper + "/environment.yaml"):
        copy(wrapper, "environment.yaml")
    print(f"Copied {wrapper} into {dst}")


# The main function to run the tests
def run(wrapper, cmd, check_log=None, extra_wrappers=[]):

    origdir = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        dst = os.path.join(d, "main")
        os.makedirs(dst, exist_ok=True)

        copy_wrapper(wrapper, dst)
        # if a wrapper depends on other wrappers, we need to include them
        for extra_wrapper in extra_wrappers:
            copy_wrapper(extra_wrapper, dst)

        # if test did not changed, not need to be run on CI action
        if (DIFF_MASTER or DIFF_LAST_COMMIT) and not any(f.startswith(wrapper) for f in DIFF_FILES):
            raise Skipped("wrappers not modified")

        # if test modified but cannot be run on CI action, raise exception
        if wrapper in WRAPPERS_TO_SKIP_ON_GA and "GITHUB_ACTION" in os.environ:
            raise Skipped("Test not runnable on GITHUB CI Action")

        # copy wrapper test and run it
        testdir = os.path.join(d, "test")
        shutil.copytree(os.path.join(wrapper, "test"), testdir)
        os.chdir(testdir)
        if os.path.exists(".snakemake"):
            shutil.rmtree(".snakemake")
        cmd = cmd + ["--wrapper-prefix", "file://{}/".format(d), "--conda-cleanup-pkgs"]

        try:
            subprocess.check_call(cmd)
        except Exception as e:
            # go back to original directory
            os.chdir(origdir)
            logfiles = [os.path.join(d, f) for d, _, files in os.walk(os.path.join(testdir, "logs")) for f in files]
            for path in logfiles:
                with open(path) as f:
                    msg = "###### Logfile: " + path + " ######"
                    print(msg, "\n")
                    print(f.read())
                    print("#" * len(msg))
            if check_log is not None:
                for f in logfiles:
                    check_log(open(f).read())
            else:
                raise e
        finally:
            # cleanup environments to save disk space
            subprocess.check_call(
                r"for env in `conda env list | grep -P '\.snakemake/conda' | "
                "cut -f1 | tr -d ' '`; do conda env remove --prefix $env; done",
                shell=True,
            )
            # go back to original directory
            os.chdir(origdir)


@skip_if_not_modified
def test_canu():
    run("wrappers/canu", ["snakemake", "--cores", "2", "--use-conda", "-F"])


@skip_if_not_modified
def test_busco():
    run(
        "wrappers/busco",
        [
            "snakemake",
            "--cores",
            "1",
            "txome_busco",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_fastqc():
    run(
        "wrappers/fastqc",
        [
            "snakemake",
            "--cores",
            "2",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_fastp_pe():
    run(
        "wrappers/fastp",
        [
            "snakemake",
            "--cores",
            "1",
            "trimmed/pe/a.1.fastq",
            "trimmed/pe/a.2.fastq",
            "report/pe/a.html",
            "report/pe/a.json",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_fastp_pe_wo_trimming():
    run(
        "wrappers/fastp",
        [
            "snakemake",
            "--cores",
            "1",
            "report/pe_wo_trimming/a.html",
            "report/pe_wo_trimming/a.json",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_fastp_se():
    run(
        "wrappers/fastp",
        [
            "snakemake",
            "--cores",
            "1",
            "trimmed/se/a.fastq",
            "report/se/a.html",
            "report/se/a.json",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_bowtie2_align():
    run(
        "wrappers/bowtie2/align",
        ["snakemake", "--cores", "1", "mapped/a.bam", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_bowtie2_build():
    run(
        "wrappers/bowtie2/build",
        ["snakemake", "--cores", "1", "genome.1.bt2", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_blast_makeblastdb_nucleotide():
    run(
        "wrappers/blast/makeblastdb",
        [
            "snakemake",
            "--cores",
            "1",
            "results/genome.fasta.ndb",
            "results/genome.fasta.nhr",
            "results/genome.fasta.nin",
            "results/genome.fasta.not",
            "results/genome.fasta.nsq",
            "results/genome.fasta.ntf",
            "results/genome.fasta.nto",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_blast_makeblastdb_protein():
    run(
        "wrappers/blast/makeblastdb",
        [
            "snakemake",
            "--cores",
            "1",
            "results/protein.fasta.pdb",
            "results/protein.fasta.phr",
            "results/protein.fasta.pin",
            "results/protein.fasta.pot",
            "results/protein.fasta.psq",
            "results/protein.fasta.ptf",
            "results/protein.fasta.pto",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_multiqc():
    run(
        "wrappers/multiqc",
        ["snakemake", "--cores", "1", "qc/multiqc.html", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_transdecoder_longorfs():
    run(
        "wrappers/transdecoder/longorfs",
        [
            "snakemake",
            "--cores",
            "1",
            "test.fa.transdecoder_dir/longest_orfs.pep",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_transdecoder_predict():
    run(
        "wrappers/transdecoder/predict",
        ["snakemake", "--cores", "1", "test.fa.transdecoder.gff3", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_trinity():
    run(
        "wrappers/trinity",
        [
            "snakemake",
            "--cores",
            "1",
            "trinity_out_dir/Trinity.fasta",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_falco():
    run(
        "wrappers/falco",
        [
            "snakemake",
            "--cores",
            "2",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_gz_to_bz2():
    run(
        "wrappers/falco",
        [
            "snakemake",
            "--cores",
            "2",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_bz2_to_gz():
    run(
        "wrappers/falco",
        [
            "snakemake",
            "--cores",
            "2",
            "--use-conda",
            "-F",
        ],
    )


@skip_if_not_modified
def test_gz_to_dsrc():
    run(
        "wrappers/falco",
        [
            "snakemake",
            "--cores",
            "2",
            "--use-conda",
            "-F",
        ],
    )


# this wrapper test do not use conda but damona
@skip_if_not_modified
def test_bcl2fastq():
    run(
        "wrappers/bcl2fastq",
        [
            "snakemake",
            "--cores",
            "2",
            "-F",
        ],
    )


@skip_if_not_modified
def test_bowtie1_build():
    run(
        "wrappers/bowtie1/build",
        ["snakemake", "--cores", "1", "--use-conda", "-F"],
    )


@skip_if_not_modified
def test_bowtie1_align():
    run(
        "wrappers/bowtie1/align",
        ["snakemake", "--cores", "1", "--use-conda", "-F"],
        extra_wrappers=["wrappers/bowtie1/build"],
    )


@skip_if_not_modified
def test_bwa_build():
    run(
        "wrappers/bwa/build",
        command,
    )


@skip_if_not_modified
def test_bwa_align():
    run("wrappers/bwa/align", command, extra_wrappers=["wrappers/bwa/build"])


@skip_if_not_modified
def test_minimap2():
    run(
        "wrappers/minimap2",
        command,
    )


@skip_if_not_modified
def test_fastq_stats():
    run(
        "wrappers/fastq_stats",
        command,
    )


@skip_if_not_modified
def test_bamtools_sort():
    run(
        "wrappers/bamtools/sort",
        command,
    )


@skip_if_not_modified
def test_bamtools_index():
    run(
        "wrappers/bamtools/index",
        command,
    )


@skip_if_not_modified
def test_feature_counts():
    run("wrappers/feature_counts", command)


@skip_if_not_modified
def test_add_read_group():
    run("wrappers/add_read_group", command)


@skip_if_not_modified
def test_snpeff_add_locus_in_fasta():
    run("wrappers/snpeff_add_locus_in_fasta", command)


@skip_if_not_modified
def test_sambamba_markdup():
    run("wrappers/sambamba_markdup", command)


@skip_if_not_modified
def test_sambamba_filter():
    run("wrappers/sambamba_filter", command)


@skip_if_not_modified
def test_samtools_depth():
    run("wrappers/samtools_depth", command)


@skip_if_not_modified
def test_sequana_coverage():
    run("wrappers/sequana_coverage", command)


@skip_if_not_modified
def test_freebayes():
    run("wrappers/freebayes", command)


@skip_if_not_modified
def test_freebayes_vcf_filter():
    run("wrappers/freebayes_vcf_filter", command)


@skip_if_not_modified
def test_snpeff():
    run("wrappers/snpeff", command)


@skip_if_not_modified
def test_rulegraph():
    run("wrappers/rulegraph", command)


@skip_if_not_modified
def test_digital_normalisation():
    run("wrappers/digital_normalisation", command)
