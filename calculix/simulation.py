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
        mesh_name="Solid-Nodes-Mesh", interface_name="Solid-Interface",
        read_data="DisplacementDelta", write_data="Force"):
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


def cleanup(sim_folder, remove_spooles=True):
    """Method to clean up all calculix files except paraview output."""
    sim_folder = Path(sim_folder)
    spooles_file = Path(__file__).resolve().parent / "spooles.out"

    if remove_spooles and spooles_file.exists():
        spooles_file.unlink()

    if sim_folder.exists() and sim_folder.is_dir():
        shutil.rmtree(sim_folder)
