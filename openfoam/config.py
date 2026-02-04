"""Module to read and parse precice configuration XML files."""
from xml.etree import ElementTree as ET


def read_precice_config(precice_cfg_path):  # noqa: WPS210
    """Reads and parses a preCICE configuration XML file.

    Args:
        precice_cfg_path (str or Path): Path to the preCICE configuration file (.xml).

    Returns:
        dict: Dictionary containing parsed configuration values.
    """
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

    ns_defs = " ".join([
        f'xmlns:{key}="{key_value}"'
        for key, key_value
        in namespaces.items()
    ])
    if "<precice-configuration " in xml_text:
        xml_text = xml_text.replace(
            "<precice-configuration", f"<precice-configuration {ns_defs}"
        )
    else:
        xml_text = xml_text.replace(
            "<precice-configuration>", f"<precice-configuration {ns_defs}>"
        )

    root = ET.fromstring(xml_text)
    dt_element = root.find(".//{*}time-window-size")

    return {"dt": dt_element.get("value")}


def read_config(filename):
    """Parses fluid config file."""
    config = {}
    with open(filename, "r", encoding="utf-8") as config_file:
        for line in config_file:
            line = line.split('#')[0].strip()
            if not line or "=" not in line:
                continue

            _add_key_value_pair(config, line)
    return config


def _add_key_value_pair(config, line):
    key, key_value = line.split('=', 1)
    key, key_value = key.strip(), key_value.strip()

    if key_value.lower() == "true":
        config[key] = True
    elif key_value.lower() == "false":
        config[key] = False
    else:
        try:
            config[key] = (
                int(float(key_value))
                if float(key_value).is_integer()
                else float(key_value)
            )
        except ValueError:
            config[key] = key_value
