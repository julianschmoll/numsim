"""Module to run CalculiX simulations with preCICE coupling and manage files."""
from pathlib import Path
import subprocess
from contextlib import contextmanager
import shutil
import yaml


@contextmanager
def precice_env_setup(environment):
    """Creates the config files required by the CalculiX adapter.

    Args:
        environment (dict): Dictionary containing the environment variables
            required to set up the preCICE configuration.

    Yields:
        None
    """
    work_dir = Path(environment.get("work_dir", Path.cwd()))
    yaml_target = work_dir / "config.yml"
    source_xml = Path(environment["cfg_path"]).resolve()
    yaml_data = {
        "participants": {
            environment["participant"]: {
                "interfaces": [{
                    "nodes-mesh": environment["mesh_name"],
                    "patch": environment["interface"],
                    "read-data": [environment["read_data"]],
                    "write-data": [environment["write_data"]],
                }]
            }
        },
        "precice-config-file": str(source_xml)
    }

    work_dir.mkdir(parents=True, exist_ok=True)
    with open(yaml_target, "w", encoding="utf-8") as yaml_file:
        yaml.dump(yaml_data, yaml_file, sort_keys=False, default_flow_style=False)
    try:
        yield
    finally:
        if yaml_target.exists():
            yaml_target.unlink()


def run(inp_file, precice_cfg, participant="Solid", ccx_cmd="ccx_preCICE"):
    """Runs a CalculiX simulation with preCICE coupling.

    Args:
        inp_file (str or Path): Path to the CalculiX input file (.inp).
        precice_cfg (str or Path): Path to the preCICE configuration file (.xml).
        participant (str): Name of the preCICE participant.
        ccx_cmd (str): Command to run CalculiX with preCICE support.

    Returns:
        Output of CalculiX simulation (.frd file path).
    """
    inp_path = Path(inp_file).resolve()
    work_dir = inp_path.parent
    inp_stem_name = inp_path.with_suffix("").name

    environment = {
        "cfg_path": precice_cfg,
        "participant": participant,
        "mesh_name": "Solid-Mesh",
        "interface": "Solid-Interface",
        "read_data": "Force",
        "write_data": "DisplacementDelta",
        "work_dir": work_dir,
    }

    with precice_env_setup(environment):
        run_cmd = [
            ccx_cmd,
            "-i", inp_stem_name,
            "-precice-participant", participant
        ]
        subprocess.run(run_cmd, check=True, cwd=work_dir)

    return work_dir / f"{inp_stem_name}.frd"


def cleanup(sim_folder, exclude=None):
    """Method to clean up all calculix files except specifically excluded.

    Args:
        sim_folder (str or Path): Path to the simulation folder.
        exclude (str or list, optional): File or list of files to exclude from deletion.

    """
    sim_folder = Path(sim_folder)

    patterns = [
        "precice*",
        "exchange*",
        "m2n*",
        "*.nam",
        "*.sur",
        "*.log",
        "*.lock",
        "spooles.out",
        "*.cvg",
        "*.dat",
        "*.inp",
        "*.sta",
        "*.12d"
    ]

    if not exclude:
        exclude = []

    if isinstance(exclude, str):
        exclude = [exclude]

    if not isinstance(exclude, list):
        raise TypeError("exclude must be a string or list")

    for pattern in patterns:
        _unlink_files(exclude, pattern, sim_folder)


def _unlink_files(exclude, pattern, sim_folder):
    """Helper to unlink files matching a pattern except those excluded.

    Args:
        exclude (list): List of filenames to exclude from deletion.
        pattern (str): Glob pattern to match files.
        sim_folder (Path): Path to the simulation folder.
    """
    for found_files in sim_folder.glob(pattern):
        if found_files.name in exclude:
            continue
        if found_files.is_dir():
            shutil.rmtree(found_files, ignore_errors=True)
        else:
            found_files.unlink()
