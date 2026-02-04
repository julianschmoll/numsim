"""Entry point for solid simulation with OpenFOAM and preCICE."""
import argparse
from pathlib import Path
import logging

import case
import config


DEFAULT_DT = 0.001


def main(scenario_cfg, precice_cfg_path):
    """Entry point for solid simulation with OpenFOAM and preCICE.

    Args:
        scenario_cfg (str or Path): Path to the scenario configuration file.
        precice_cfg_path (str or Path): Path to the preCICE configuration file.
    """
    simulation_folder = Path(__file__).resolve().parent / "out"

    cfg = config.read_config(scenario_cfg)
    cfg["coupled"] = False

    if precice_cfg_path:
        cfg["coupled"] = True
        cfg["precice_cfg"] = str(Path(precice_cfg_path).resolve())
        dt = config.read_precice_config(precice_cfg_path).get("dt", DEFAULT_DT)
        cfg["maximumDt"] = float(dt)

    case.generate(simulation_folder, cfg)
    case.run(
        str(simulation_folder.resolve()),
        Path(scenario_cfg).stem
    )


# entry point for solid simulation with openfoam and precice
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
        default=(Path(__file__).resolve().parent.parent
                 / "cfg" / "fluid" / "lid_driven_cavity.txt"),
        help="Scenario to run with calculix",
    )

    parser.add_argument(
        "--cleanup",
        type=bool,
        default=True,
        help="Clean up simulation files after run",
    )

    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    main(args.scenario, args.precice_cfg)
