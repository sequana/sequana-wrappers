import importlib


def get_shell(tool_path: str, version: str) -> str:
    """Return a shell command string from the sequana_wrappers shell library.

    ``tool_path`` is a slash-separated string encoding the tool and command
    (and any future sub-commands), e.g. ``"bwa/align"`` or ``"bamtools/stats"``.

    Every version — including ``"dev"`` — maps to a subdirectory of the same
    name under ``shells/<tool>/<command>/``.  No silent fallback: an explicit
    error is raised if the requested version does not exist.

    :param tool_path: slash-separated tool/command path, e.g. ``"bwa/align"``
    :param version: shell command version, e.g. ``"v1"``, or ``"dev"`` for the
        development (unreleased) version.

    Example usage in a pipeline rules file::

        shell: manager.get_shell("bwa/align", "v1")    # pinned — recommended
        shell: manager.get_shell("bwa/align", "dev")   # development version
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
