"""Command-line interface for sequana_wrappers."""

import argparse
import subprocess
from pathlib import Path


def get_version() -> str:
    """Get package version from pyproject.toml."""
    pyproject = Path(__file__).parent.parent / "pyproject.toml"
    with open(pyproject) as f:
        for line in f:
            if line.startswith('version = "'):
                return line.split('"')[1]
    return "unknown"


def count_git_tags() -> int:
    """Count git tags."""
    try:
        result = subprocess.run(
            ["git", "tag"],
            capture_output=True,
            text=True,
            cwd=Path(__file__).parent.parent,
        )
        return len(result.stdout.strip().split("\n")) if result.stdout.strip() else 0
    except Exception:
        return 0


def count_shells() -> int:
    """Count shells in sequana_wrappers/shells/."""
    shells_dir = Path(__file__).parent / "shells"
    return len([d for d in shells_dir.iterdir() if d.is_dir()])


def count_wrappers() -> int:
    """Count wrappers in wrappers/."""
    wrappers_dir = Path(__file__).parent.parent / "wrappers"
    return len([d for d in wrappers_dir.iterdir() if d.is_dir()])


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="sequana_wrappers command-line interface",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"sequana_wrappers {get_version()}",
    )
    parser.add_argument(
        "--stats",
        action="store_true",
        help="Show library statistics",
    )

    args = parser.parse_args()

    if args.stats:
        version = get_version()
        tags = count_git_tags()
        shells = count_shells()
        wrappers = count_wrappers()

        print(f"sequana_wrappers {version}")
        print(f"  Git tags:      {tags}")
        print(f"  Shells:        {shells}")
        print(f"  Wrappers:      {wrappers} (legacy, deprecated)")
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
