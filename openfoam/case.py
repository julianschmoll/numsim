from pathlib import Path

import subprocess
import logging
from string import Template

import math


def get_template(name, coupled=False):
    template_folder = (
            Path(__file__).parent.parent / "resources" / "openfoam"
    )
    template_path = template_folder / name

    if coupled:
        precice_template = template_folder / f"{name}Precice"
        template_path = precice_template \
            if precice_template.exists() \
            else template_path

    with open(template_path, "r") as template_file:
        return Template(template_file.read())


def add_openfoam_keys(cfg):
    """Calculates physical properties and configures boundary conditions."""
    walls = ["Bottom", "Top", "Left", "Right"]
    has_outflow_boundary = False
    u_max = 0.0

    for wall in walls:
        b_type = cfg.get(f"boundary{wall}", "")
        vx = cfg.get(f"dirichlet{wall}X", 0.0)
        vy = cfg.get(f"dirichlet{wall}Y", 0.0)
        u_max = max(u_max, (vx**2 + vy**2)**0.5)

        point_cond = "type fixedValue; value uniform (0 0 0);"

        if "Coupled" in b_type and cfg.get("coupled"):
            u_cond = "type movingWallVelocity; value uniform (0 0 0);"
            p_cond = "type zeroGradient;"
            point_cond = "type fixedValue; value $internalField;"
        elif "Outflow" in b_type:
            u_cond = "type zeroGradient;"
            p_cond = "type fixedValue; value uniform 0;"
            has_outflow_boundary = True
        elif vx or vy:
            if "frequency" in cfg:
                cfg["scale"] = 2 * math.pi * cfg["frequency"] * cfg["timeShift"]
                cfg["level"] = f"({vx} {vy} {0})"
                u_cond = get_template("dynamicBoundary").substitute(cfg)
                p_cond = "type fixedFluxPressure; value uniform 0;"
            else:
                u_cond = f"type fixedValue; value uniform ({vx} {vy} 0);"
                p_cond = "type zeroGradient;"
        else:
            u_cond = "type noSlip;"
            p_cond = "type zeroGradient;"

        cfg[f"u{wall}"] = u_cond
        cfg[f"p{wall}"] = p_cond
        cfg[f"point{wall}"] = point_cond

    # ToDo: This is uggo
    u_ref = u_max if u_max > 0 else 1.0
    l_ref = cfg.get('characteristicLength', cfg.get('physicalSizeY', 2.0))
    cfg["nu"] = (u_ref * l_ref) / cfg['re']

    if has_outflow_boundary:
        cfg["pRefEntry"] = "// pRef driven by Outflow boundary"
    else:
        cfg["pRefEntry"] = "pRefCell 0; pRefValue 0;"


def add_precice_keys(cfg):
    walls = ["Bottom", "Top", "Left", "Right"]
    cfg["coupled_patches"] = " ".join(
        wall.lower() for wall in walls
        if "Coupled" in cfg.get(f"boundary{wall}", "")
    )


def write_file(path_str, content):
    path = Path(path_str)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def generate(simulation_folder, cfg):
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
    cmd_str = " ".join(str(substr) for substr in cmd)
    logging.info(f"Running: {cmd_str}")
    subprocess.run(cmd, check=True, cwd=case_dir)


def run(case_path, case_name="test"):
    _run_command(["blockMesh", "-case", case_path], case_path)
    _run_command(["checkMesh", "-case", case_path], case_path)
    _run_command(["pimpleFoam", "-case", case_path], case_path)
    _run_command(["touch", f"{case_name}.foam", case_path], case_path)
