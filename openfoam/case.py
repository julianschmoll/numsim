from pathlib import Path

import subprocess
import logging
from string import Template


def get_template(name):
    """Helper to get template content."""
    template_path = (
            Path(__file__).parent.parent / "resources" / "openfoam" / name
    )
    with open(template_path, "r") as template_file:
        return template_file.read()


def write_file(path_str, content):
    """Helper to write content to a file, creating parent directories if needed."""
    path = Path(path_str)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def generate_null(simulation_folder, cfg):
    write_file(
        simulation_folder / "0" / "U",
        get_template("U")
    )
    write_file(
        simulation_folder / "0" / "p",
        get_template("p")
    )


def generate_constant(simulation_folder, cfg):

    constant_folder = Path(simulation_folder) / "constant"
    constant_files = ["transportProperties", "turbulenceProperties"]

    for filename in constant_files:
        write_file(
            constant_folder / filename,
            Template(get_template(filename)).substitute(cfg)
        )


def generate_system(simulation_folder, cfg):
    system_folder = Path(simulation_folder) / "system"
    system_files = [
        "controlDict", "fvSchemes", "fvSolution",
        "decomposeParDict", "PDRblockMeshDict", "blockMeshDict"
    ]
    for filename in system_files:
        print(filename)
        write_file(
            system_folder / filename,
            Template(get_template(filename)).substitute(cfg)
        )


def generate(simulation_folder, cfg):
    """Generates a simple OpenFOAM case structure."""
    (simulation_folder / "0").mkdir(parents=True, exist_ok=True)
    (simulation_folder / "constant").mkdir(parents=True, exist_ok=True)
    (simulation_folder / "system").mkdir(parents=True, exist_ok=True)

    generate_null(simulation_folder, cfg)
    generate_constant(simulation_folder, cfg)
    generate_system(simulation_folder, cfg)


def _run_command(cmd, case_dir: Path):
    cmd_str = " ".join(str(substr) for substr in cmd)
    logging.info(f"Running: {cmd_str}")
    subprocess.run(cmd, check=True, cwd=case_dir)

def run(case_path):
    _run_command(["blockMesh", "-case", case_path], case_path)
    _run_command(["checkMesh", "-case", case_path], case_path)
    _run_command(["pimpleFoam", "-case", case_path], case_path)
    _run_command(["foamToVTK", "-case", case_path], case_path)
