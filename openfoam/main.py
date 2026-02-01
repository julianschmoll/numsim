import argparse
from pathlib import Path
import logging
import xml.etree.ElementTree as ET

import case

def read_precice_config(precice_cfg_path):
    namespaces = {
        'data': 'data',
        'm2n': 'm2n',
        'coupling-scheme': 'coupling-scheme',
        'participant': 'participant',
        'mesh': 'mesh',
        'mapping': 'mapping',
        'export': 'export',
        'acceleration': 'acceleration'
    }

    for prefix, uri in namespaces.items():
        ET.register_namespace(prefix, uri)

    with open(precice_cfg_path, "r", encoding="utf-8") as xml_file:
        xml_text = xml_file.read()

    ns_defs = " ".join([f'xmlns:{k}="{v}"' for k, v in namespaces.items()])
    if "<precice-configuration " in xml_text:
        xml_text = xml_text.replace("<precice-configuration", f"<precice-configuration {ns_defs}")
    else:
        xml_text = xml_text.replace("<precice-configuration>", f"<precice-configuration {ns_defs}>")

    root = ET.fromstring(xml_text)
    dt_element = root.find(".//{*}time-window-size")

    if dt_element is not None:
        return {"dt": dt_element.get("value")}


def read_config(filename):
    """Parses config. If a key appears multiple times with tuple values, it collects them."""
    config = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.split('#')[0].strip()
            if not line or "=" not in line:
                continue

            key, value = line.split('=', 1)
            key, value = key.strip(), value.strip()

            if value.lower() == "true":
                config[key] = True
            elif value.lower() == "false":
                config[key] = False
            else:
                try:
                    num_val = float(value)
                    config[key] = int(num_val) if num_val.is_integer() else num_val
                except ValueError:
                    config[key] = value
    return config

def main(scenario_cfg, precice_cfg_path, cleanup=True):
    simulation_folder = Path(__file__).resolve().parent / "out"

    cfg = read_config(scenario_cfg)
    cfg["coupled"] = False

    if precice_cfg_path:
        cfg["coupled"] = True
        cfg["precice_cfg"] = str(Path(precice_cfg_path).resolve())
        dt = read_precice_config(precice_cfg_path).get("dt", 0.001)
        cfg["maximumDt"] = float(dt)

    case.generate(simulation_folder, cfg)
    case.run(
        str(simulation_folder.resolve()),
        Path(scenario_cfg).stem
    )

    # convert to vtk
    if cleanup:
        pass


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
        default=Path(__file__).resolve().parent.parent
                 / "cfg" / "fluid" / "lid_driven_cavity.txt",
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

    main(args.scenario, args.precice_cfg, cleanup=args.cleanup)
