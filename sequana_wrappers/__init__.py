import importlib


def get_shell(tool_path: str, version: str) -> str:
    """Return a shell command string from the sequana_wrappers shell library.

    :param tool_path: slash-separated tool/command path, e.g. ``"bwa/align"``
    :param version: shell command version, e.g. ``"v1"``, or ``"dev"`` for the
        development (unreleased) version.

    Example usage in a pipeline rules file::

        shell: manager.get_shell("bwa/align", "v1")
    """
    parts = tool_path.split("/")
    module_path = ".".join(["sequana_wrappers", "shells"] + parts + [version, "cmd"])
    try:
        return importlib.import_module(module_path).CMD
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            f"Shell command version '{version}' not found: '{module_path}'.\n"
            f"Available versions are the subdirectories under "
            f"sequana_wrappers/shells/{'/'.join(parts)}/.\n"
            f"Use version='dev' for the development version."
        )


def get_run(tool_path: str, version: str):
    """Return a Python callable from the sequana_wrappers snippets library.

    The returned callable has the signature ``execute(snakemake)`` and is
    intended for use inside Snakemake ``run:`` blocks.

    :param tool_path: slash-separated tool/command path, e.g. ``"rulegraph/run"``
    :param version: snippet version, e.g. ``"v1"``, or ``"dev"`` for the
        development (unreleased) version.

    Example usage in a pipeline rules file::

        run:
            manager.get_run("rulegraph/run", "v1")(snakemake)
    """
    parts = tool_path.split("/")
    module_path = ".".join(["sequana_wrappers", "snippets"] + parts + [version, "code"])
    try:
        return importlib.import_module(module_path).execute
    except ModuleNotFoundError:
        raise ModuleNotFoundError(
            f"Snippet version '{version}' not found: '{module_path}'.\n"
            f"Available versions are the subdirectories under "
            f"sequana_wrappers/snippets/{'/'.join(parts)}/.\n"
            f"Use version='dev' for the development version."
        )
