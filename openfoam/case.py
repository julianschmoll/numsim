"""Generates OpenFOAM case files based on configuration."""
from pathlib import Path

import subprocess
import logging
from string import Template

import math

ZERO_GRADIENT = "type zeroGradient;"


def get_template(name, coupled=False):
    """Retrieves the template file for OpenFOAM case generation.

    Args:
        name (str): Name of the template file.
        coupled (bool): Whether to use the preCICE-coupled template.

    Returns:
        Template: The template object for the specified file.
    """
    template_folder = (
            Path(__file__).parent.parent / "resources" / "openfoam"
    )
    template_path = template_folder / name

    if coupled:
        precice_template = template_folder / f"{name}Precice"
        template_path = precice_template \
            if precice_template.exists() \
            else template_path

    with open(
            template_path, "r", encoding="utf-8"
    ) as template_file:
        return Template(template_file.read())


def add_openfoam_keys(cfg):  # noqa: WPS210, WPS231
    """Calculates physical properties and configures boundary conditions.

    Args:
        cfg (dict): Configuration dictionary to be updated.
    """
    walls = ["Bottom", "Top", "Left", "Right"]
    has_outflow_boundary = False
    u_max = .0

    for wall in walls:
        b_type = cfg.get(f"boundary{wall}", "")
        vx = cfg.get(f"dirichlet{wall}X", .0)
        vy = cfg.get(f"dirichlet{wall}Y", .0)
        u_max = max(u_max, (vx**2 + vy**2)**0.5)

        point_cond = "type fixedValue; value uniform (0 0 0);"

        if "Coupled" in b_type and cfg.get("coupled"):
            u_cond = "type movingWallVelocity; value uniform (0 0 0);"
            p_cond = ZERO_GRADIENT
            point_cond = "type fixedValue; value $internalField;"
        elif "Outflow" in b_type:
            u_cond = ZERO_GRADIENT
            p_cond = "type fixedValue; value uniform 0;"
            has_outflow_boundary = True
        elif vx or vy:
            if "startBurst" in cfg and "endBurst" in cfg:
                cfg["initialU"] = f"({vx} {vy} 0)" \
                    if cfg["startBurst"] <= 0 \
                    else "(0 0 0)"
                cfg["u"] = f"({vx} {vy} 0)"
                u_cond = get_template("burstBoundary").substitute(cfg)
                p_cond = "type fixedFluxPressure; value uniform 0;"
            elif "frequency" in cfg:
                cfg["scale"] = (2 * math.pi
                                * cfg.get("frequency") * cfg.get("timeShift"))
                cfg["level"] = f"({vx} {vy} {0})"  # noqa: WPS237
                u_cond = get_template("dynamicBoundary").substitute(cfg)
                p_cond = "type fixedFluxPressure; value uniform 0;"
            else:
                u_cond = f"type fixedValue; value uniform ({vx} {vy} 0);"
                p_cond = ZERO_GRADIENT
        else:
            u_cond = "type noSlip;"
            p_cond = ZERO_GRADIENT

        cfg[f"u{wall}"] = u_cond
        cfg[f"p{wall}"] = p_cond
        cfg[f"point{wall}"] = point_cond

    u_ref = u_max if u_max > 0 else 1.0
    l_ref = cfg.get('characteristicLength', cfg.get('physicalSizeX', 2))
    cfg["nu"] = (u_ref * l_ref) / cfg['re']

    if has_outflow_boundary:
        cfg["pRefEntry"] = "// pRef driven by Outflow boundary"
    else:
        cfg["pRefEntry"] = "pRefCell 0; pRefValue 0;"


def add_precice_keys(cfg):
    """Adds preCICE-related keys to the configuration.

    Args:
        cfg (dict): Configuration dictionary to be updated.
    """
    walls = ["Bottom", "Top", "Left", "Right"]
    cfg["coupled_patches"] = " ".join(
        wall.lower() for wall in walls
        if "Coupled" in cfg.get(f"boundary{wall}", "")
    )


def write_file(path_str, file_content):
    """Writes content to a file, creating parent directories if needed.

    Args:
        path_str (str or Path): Path to the file.
        file_content (str): Content to write to the file.
    """
    path = Path(path_str)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(file_content, encoding="utf-8")


def generate(simulation_folder, cfg):
    """Generates OpenFOAM case files based on configuration.

    Args:
        simulation_folder (str or Path): Path to the simulation folder.
        cfg (dict): Configuration dictionary.
    """
    case_structure = {
        "0": ["U", "p"],
        "constant": ["transportProperties", "turbulenceProperties"],
        "system": [
            "controlDict", "fvSchemes", "fvSolution",
            "decomposeParDict", "PDRblockMeshDict", "blockMeshDict"
        ]
    }

    add_openfoam_keys(cfg)

    if cfg.get("coupled"):
        case_structure["0"].append("pointDisplacement")
        case_structure["constant"].append("dynamicMeshDict")
        case_structure["system"].append("preciceDict")
        add_precice_keys(cfg)

    for folder, files in case_structure.items():
        for filename in files:
            write_file(
                simulation_folder / folder / filename,
                get_template(filename, cfg.get("coupled")).substitute(cfg)
            )


def _run_command(cmd, case_dir: Path):
    """Runs a command in a specified directory.

    Args:
        cmd (list): Command and its arguments as a list.
        case_dir (str or Path): Directory to run the command in.
    """
    cmd_str = " ".join(str(substr) for substr in cmd)
    logging.info(f"Running: {cmd_str}")
    subprocess.run(cmd, check=True, cwd=case_dir)


def run(case_path, case_name="test"):
    """Runs the OpenFOAM simulation.

    Args:
        case_path (str or Path): Path to the OpenFOAM case directory.
        case_name (str): Name for the output .foam file.
    """
    _run_command(["blockMesh", "-case", case_path], case_path)
    _run_command(["checkMesh", "-case", case_path], case_path)
    _run_command(["pimpleFoam", "-case", case_path], case_path)
    _run_command(["touch", f"{case_name}.foam", case_path], case_path)
