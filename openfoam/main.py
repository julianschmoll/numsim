import argparse
import logging
import subprocess
from pathlib import Path

import generate, convert


def read_config(filename):
    """Parses the custom key=value config file."""
    config = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line:
                continue
            if '=' in line:
                key, value = line.split('=', 1)
                key, value = key.strip(), value.strip()
                if value.lower() == 'true':
                    config[key] = True
                elif value.lower() == 'false':
                    config[key] = False
                else:
                    try:
                        config[key] = float(value)
                        if config[key].is_integer():
                            config[key] = int(config[key])
                    except ValueError:
                        config[key] = value
    return config


def generate_case(cfg_path: Path, case_dir: Path):
    config = read_config(cfg_path)
    generate.write_file(
        case_dir / "system" / "blockMeshDict",
        generate.generate_block_mesh(config)
    )
    generate.write_file(
        case_dir / "system" / "controlDict",
        generate.generate_control_dict(config)
    )
    generate.write_file(
        case_dir / "system" / "fvSchemes",
        generate.generate_fv_schemes(config)
    )
    generate.write_file(
        case_dir / "system" / "fvSolution",
        generate.generate_fv_solution(config)
    )
    generate.write_file(
        case_dir / "system" / "decomposeParDict",
        generate.generate_decompose_par_dict(config)
    )
    generate.write_file(
        case_dir / "system" / "preciceDict",
        generate.generate_precice_dict(config)
    )

    # constant
    generate.write_file(
        case_dir / "constant" / "transportProperties",
        generate.generate_transport_properties(config)
    )
    generate.write_file(
        case_dir / "constant" / "turbulenceProperties",
        generate.generate_turbulence_properties(config)
    )
    generate.write_file(
        case_dir / "constant" / "dynamicMeshDict",
        generate.generate_dynamic_mesh_dict(config)
    )

    # fields
    generate.write_file(
        case_dir / "0" / "U",
        generate.generate_field(config, "U")
    )
    generate.write_file(
        case_dir / "0" / "p",
        generate.generate_field(config, "p")
    )
    generate.write_file(
        case_dir / "0" / "pointDisplacement",
        generate.generate_field(config, "pointDisplacement")
    )


def _run_command(cmd, case_dir: Path):
    cmd_str = " ".join(str(substr) for substr in cmd)
    logging.info(f"Running: {cmd_str}")
    subprocess.run(cmd, check=True, cwd=case_dir)


def main(config_file: str, case_dir: str, cleanup: bool = True):
    cfg_path = Path(config_file).resolve()
    case_path = Path(case_dir).resolve()
    case_path.mkdir(parents=True, exist_ok=True)

    logging.info("Generating OpenFOAM case in %s using config %s", case_path, cfg_path)
    generate_case(cfg_path, case_path)

    _run_command(["blockMesh", "-case", case_path], case_path)
    _run_command(["checkMesh", "-case", case_path], case_path)
    _run_command(["pimpleFoam", "-case", case_path], case_path)

    out_dir = convert.export_vtk(case_path)
    logging.info("VTK files exported to %s", out_dir)

    if cleanup:
        raise NotImplementedError("Cleanup not implemented for OpenFOAM cases yet.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser("OpenFOAM fluid participant entrypoint.")
    parser.add_argument(
        "--config",
        type=str,
        default=None,
        help="Path to case config file."
    )
    parser.add_argument(
        "--case",
        type=str,
        default=Path(__file__).resolve().parent / "case",
        help="Directory where case files and run will be placed."
    )
    parser.add_argument(
        "--cleanup",
        type=bool,
        default=False,
        help="Whether to perform cleanup after run."
    )

    args = parser.parse_args()
    cfg = args.config or (Path(__file__).resolve().parent / "lid_driven_cavity.txt")
    logging.basicConfig(level=logging.INFO)
    main(cfg, args.case, cleanup=args.cleanup)
