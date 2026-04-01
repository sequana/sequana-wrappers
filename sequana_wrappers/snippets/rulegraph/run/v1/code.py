import os
from pathlib import Path

from sequana_pipetools import SequanaConfig
from sequana_pipetools.snaketools import DOTParser


def execute(input, output, params):
    input_filename = input[0]
    output_dot = output[0]
    config_filename = params.get("configname")
    mapper_urls = params.get("mapper", {})
    required_local_files = params.get("required_local_files", [])

    def parse_path(dico):
        for key, value in dico.items():
            try:
                if value:
                    try:
                        dico[key] = str(Path(value).resolve(strict=True))
                    except FileNotFoundError:
                        pass
            except (TypeError, OverflowError):
                try:
                    parse_path(value)
                except AttributeError:
                    pass

    cfg = SequanaConfig(config_filename)
    parse_path(cfg.config)
    cfg._update_yaml()

    def link_required_files(required_paths):
        for path in required_paths:
            try:
                if Path(path).exists():
                    Path(f"rulegraph/{path}").symlink_to(f"../{path}")
                else:
                    raise FileNotFoundError(f"Required local file {path} is not found.")
            except FileExistsError:
                pass

    cwd = Path.cwd()
    try:
        Path("rulegraph").mkdir(exist_ok=True)
        if config_filename:
            cfg.save(filename=Path("rulegraph") / config_filename)
        link_required_files(required_local_files)
        os.system(f'cd rulegraph && snakemake -s "{input_filename}" --rulegraph --nolock > rg.dot; cd ..')
    except Exception as err:
        print(err)
        os.chdir(cwd)

    d = DOTParser(cwd / "rulegraph" / "rg.dot")
    d.add_urls(output_filename=output_dot, mapper=mapper_urls)
