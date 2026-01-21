import subprocess

from pathlib import Path
import shutil


def run(inp_file, precice_cfg, cxx_cmd="ccx"):
    run_cmd = [
        cxx_cmd,
        inp_file,
    ]
    subprocess.Popen(run_cmd, stdout=subprocess.PIPE).communicate()

    return f"{inp_file}.frd"


def cleanup(sim_folder, remove_spooles=True):
    """Method to clean up all calculix files except paraview output."""
    sim_folder = Path(sim_folder)
    spooles_file = Path(__file__).resolve().parent / "spooles.out"

    if remove_spooles and spooles_file.exists():
        spooles_file.unlink()

    if sim_folder.exists():
        shutil.rmtree(sim_folder)
