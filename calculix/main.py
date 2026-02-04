"""Entrypoint for solid simulation with calculix and precice."""
import argparse
from pathlib import Path
import logging

from geometry import Geometry
import simulation
import writer


def main(scenario_cfg, precice_cfg_path, cleanup=True):
    """Entrypoint for solid simulation with calculix and precice.

    Args:
        scenario_cfg (str or Path): Path to the scenario configuration file.
        precice_cfg_path (str or Path): Path to the preCICE configuration file.
        cleanup (bool, optional): Whether to clean up simulation files after run.
    """
    simulation_folder = Path(__file__).resolve().parent / "out"
    geometry = Geometry(scenario_cfg)

    inp_path = geometry.write_file(
        simulation_folder / "geo.inp",
        mesh_name="Solid-Mesh",
        interface_name="Solid-Interface",
        )

    sim_out = simulation.run(
        inp_path,
        precice_cfg_path,
    )
    writer.convert_to_vtk(sim_out, "out/solid.vtk")
    if cleanup:
        simulation.cleanup(simulation_folder)


# entry point for solid simulation with calculix and precice
if __name__ == "__main__":
    parser = argparse.ArgumentParser("Solid participant entrypoint.")
    parser.add_argument(
        "--precice_cfg",
        type=str,
        default=None,
        help="Path to the preCICE config file"
    )
    parser.add_argument(
        "--scenario",
        type=str,
        default=None,
        help="Scenario to run with calculix",
    )

    parser.add_argument(
        "--cleanup",
        type=bool,
        default=False,
        help="Clean up simulation files after run",
    )

    args = parser.parse_args()
    precice_cfg = args.precice_cfg
    if not precice_cfg:
        precice_cfg = (Path(__file__).resolve().parent.parent
                       / "cfg" / "precice" / "2d_elastic_tube.xml")

    scenario = args.scenario
    if not scenario:
        scenario = (Path(__file__).resolve().parent.parent
                    / "cfg" / "calculix" / "2d_elastic_tube.yaml")

    logging.basicConfig(level=logging.INFO)

    main(scenario, precice_cfg, cleanup=args.cleanup)
