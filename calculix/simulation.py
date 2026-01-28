from pathlib import Path
import subprocess
from contextlib import contextmanager
import shutil
import yaml

@contextmanager
def precice_env_setup(environment):
    """Creates the config files required by the CalculiX adapter in the input file folder."""
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
    with open(yaml_target, "w") as yaml_file:
        yaml.dump(yaml_data, yaml_file, sort_keys=False, default_flow_style=False)
    try:
        yield
    finally:
        if yaml_target.exists():
            yaml_target.unlink()


def run(inp_file, precice_cfg, participant="Solid", ccx_cmd="ccx_preCICE",
        mesh_name="Solid-Mesh", interface_name="Solid-Interface",
        read_data="Force", write_data="DisplacementDelta"):
    inp_path = Path(inp_file).resolve()
    work_dir = inp_path.parent
    inp_stem_name = inp_path.with_suffix("").name

    environment = {
        "cfg_path": precice_cfg,
        "participant": participant,
        "mesh_name": mesh_name,
        "interface": interface_name,
        "read_data": read_data,
        "write_data": write_data,
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
    """Method to clean up all calculix files except specifically excluded."""
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
        for p in sim_folder.glob(pattern):
            if p.name in exclude:
                continue
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
            else:
                p.unlink()
