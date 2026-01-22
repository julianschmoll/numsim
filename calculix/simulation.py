from pathlib import Path
import subprocess
from contextlib import contextmanager
import shutil
import yaml

@contextmanager
def precice_env_setup(environment):
    """Creates the config files required by the CalculiX adapter."""
    yaml_target = Path("config.yml")
    source_xml = Path(environment["cfg_path"]).resolve()
    yaml_data = {
        "participants": {
            environment["participant"]: {
                "interfaces": [{
                    "nodes-mesh-with-connectivity": environment["mesh_name"],
                    "patch": environment["interface"],
                    "read-data": [environment["read_data"]],
                    "write-data": [environment["write_data"]],
                }]
            }
        },
        "precice-config-file": str(source_xml)
    }
    with open(yaml_target, "w") as yaml_file:
        yaml.dump(yaml_data, yaml_file, sort_keys=False, default_flow_style=False)
    try:
        yield
    finally:
        if yaml_target.exists():
            yaml_target.unlink()

def run(inp_file, precice_cfg, participant="Solid", ccx_cmd="ccx_preCICE"):
    environment = {
        "cfg_path": precice_cfg,
        "participant": participant,
        "mesh_name": "Solid-Nodes-Mesh",
        "interface": "Solid-Interface",
        "read_data": "DisplacementDelta",
        "write_data": "Force",
    }
    with precice_env_setup(environment):
        inp_stem = str(Path(inp_file).with_suffix(""))
        run_cmd = [
            ccx_cmd,
            "-i", inp_stem,
            "-precice-participant", participant
        ]
        subprocess.run(run_cmd, check=True)
    return Path(f"{inp_stem}.frd")


def cleanup(sim_folder, remove_spooles=True):
    """Method to clean up all calculix files except paraview output."""
    sim_folder = Path(sim_folder)
    spooles_file = Path(__file__).resolve().parent / "spooles.out"

    if remove_spooles and spooles_file.exists():
        spooles_file.unlink()

    if sim_folder.exists() and sim_folder.is_dir():
        shutil.rmtree(sim_folder)
