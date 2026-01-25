from pathlib import Path

import subprocess
import logging
from string import Template


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

        if "Outflow" in b_type:
            u_cond = "type zeroGradient;"
            p_cond = "type fixedValue; value uniform 0;"
            has_outflow_boundary = True
        elif vx or vy:
            u_cond = f"type fixedValue; value uniform ({vx} {vy} 0);"
            p_cond = "type zeroGradient;"
        else:
            u_cond = "type noSlip;"
            p_cond = "type zeroGradient;"

        cfg[f"u{wall}"] = u_cond
        cfg[f"p{wall}"] = p_cond

    # ToDo: This is uggo
    u_ref = u_max if u_max > 0 else 1.0
    l_ref = cfg.get('characteristicLength', cfg.get('physicalSizeY', 2.0))
    cfg["nu"] = (u_ref * l_ref) / cfg['re']

    if has_outflow_boundary:
        cfg["pRefEntry"] = "// pRef driven by Outflow boundary"
    else:
        cfg["pRefEntry"] = "pRefCell 0; pRefValue 0;"


def get_template(name):
    template_path = (
            Path(__file__).parent.parent / "resources" / "openfoam" / name
    )
    with open(template_path, "r") as template_file:
        return Template(template_file.read())


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

    for folder, files in case_structure.items():
        for filename in files:
            write_file(
                simulation_folder / folder / filename,
                get_template(filename).substitute(cfg)
            )


def _run_command(cmd, case_dir: Path):
    cmd_str = " ".join(str(substr) for substr in cmd)
    logging.info(f"Running: {cmd_str}")
    subprocess.run(cmd, check=True, cwd=case_dir)


def run(case_path):
    _run_command(["blockMesh", "-case", case_path], case_path)
    _run_command(["checkMesh", "-case", case_path], case_path)
    _run_command(["pimpleFoam", "-case", case_path], case_path)
    _run_command(["foamToVTK", "-case", case_path], case_path)
