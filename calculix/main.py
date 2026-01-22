import argparse
from pathlib import Path
import logging

from geometry import Geometry
import simulation
import writer


def main(scenario_cfg, precice_cfg_path, cleanup=True):
    simulation_folder = Path(__file__).resolve().parent / "out"
    geometry = Geometry(scenario_cfg)

    mesh_name = "Solid-Nodes-Mesh"
    interface_name = "Solid-Interface"

    inp_path = geometry.write_file(
        simulation_folder / "geo.inp",
        mesh_name=mesh_name,
        interface_name=interface_name,
        )

    sim_out = simulation.run(
        inp_path,
        precice_cfg_path,
        mesh_name=mesh_name,
        interface_name=interface_name,
    )
    writer.convert_to_vtk(sim_out, "out/output.vtk")
    if cleanup:
        simulation.cleanup(simulation_folder, remove_spooles=True)


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
        default=True,
        help="Clean up simulation files after run",
    )

    args = parser.parse_args()
    precice_cfg = args.precice_cfg
    if not precice_cfg:
        precice_cfg = (Path(__file__).resolve().parent.parent
                       / "cfg" / "precice" / "dummy_config.xml")

    scenario = args.scenario
    if not scenario:
        scenario = (Path(__file__).resolve().parent.parent
                    / "cfg" / "calculix" / "2d_elastic_tube.yaml")

    logging.basicConfig(level=logging.INFO)

    main(scenario, precice_cfg, cleanup=args.cleanup)
