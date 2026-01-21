import argparse
from pathlib import Path

import geometry
import simulation
import writer

def main(scenario_cfg, precice_cfg_path):
    geometry_file = geometry.generate(scenario_cfg)
    sim_out = simulation.run(geometry_file, precice_cfg_path)
    writer.convert_to_vtk(sim_out, "output.vtk")
    simulation.cleanup()


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
    args = parser.parse_args()

    precice_cfg = args.precice_cfg
    if not precice_cfg:
        precice_cfg = (Path(__file__).resolve().parent.parent
                            / "cfg" / "precice" / "dummy_config.xml")

    scenario = args.scenario
    if not scenario:
        scenario = (Path(__file__).resolve().parent.parent
                         / "cfg" / "calculix" / "2d_elastic_tube.yaml")

    main(scenario, precice_cfg)
