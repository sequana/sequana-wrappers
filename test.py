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
    DIFF_FILES = (
        subprocess.check_output(["git", "diff", compare, "--name-only"])
        .decode()
        .split("\n")
    )


class Skipped(Exception):
    pass


skip_if_not_modified = pytest.mark.xfail(raises=Skipped)


def run(wrapper, cmd, check_log=None):

    origdir = os.getcwd()
    with tempfile.TemporaryDirectory() as d:
        dst = os.path.join(d, "main")
        os.makedirs(dst, exist_ok=True)

        copy = lambda pth, src: shutil.copy(
            os.path.join(pth, src), os.path.join(dst, pth)
        )

        success = False
        for ext in ("py", "R", "Rmd"):
            script = "wrapper." + ext
            if os.path.exists(os.path.join(wrapper, script)):
                os.makedirs(os.path.join(dst, wrapper), exist_ok=True)
                copy(wrapper, script)
                success = True
                break
        assert success, "No wrapper script found for {}".format(wrapper)
        copy(wrapper, "environment.yaml")

        if (DIFF_MASTER or DIFF_LAST_COMMIT) and not any(f.startswith(wrapper) for f in DIFF_FILES):
            raise Skipped("wrappers not modified")

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
            logfiles = [
                os.path.join(d, f)
                for d, _, files in os.walk(os.path.join(testdir, "logs"))
                for f in files
            ]
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
    run(
        "wrappers/canu",
        [
            "snakemake",
            "--cores",
            "2",
            "--use-conda",
            "-F"
        ]
    )


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
