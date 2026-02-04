"""Generates geometry for CalculiX simulations based on config file."""
from pathlib import Path

import yaml

BOTTOM = "bottom"
TOP = "top"
WALLS = "walls"
MAX_LINE_LENGTH = 16
GEOMETRY = "geometry"
DEFAULT_DT = 0.01
DEFAULT_BETA = 0.001


def read_config(config_path):
    """Reads a YAML config and returns its contents as a dictionary.

    Args:
        config_path (str or Path): Path to the YAML configuration file.

    Returns:
        dict: Contents of the YAML configuration file.
    """
    if not Path(config_path).exists():
        raise FileNotFoundError
    if not Path(config_path).is_file():
        raise FileNotFoundError
    with open(config_path, "r", encoding="utf-8") as config_file:
        return yaml.load(config_file, Loader=yaml.FullLoader)
    raise RuntimeError("Could not read configuration file.")


def _write_node_list(open_file, nodes):
    """Helper to write a list of nodes formatted for CalculiX

    Args:
        open_file: File object to write to.
        nodes (list of int): List of node IDs.
    """
    for index in range(0, len(nodes), MAX_LINE_LENGTH):
        chunk = nodes[index:index + MAX_LINE_LENGTH]
        open_file.write(", ".join(map(str, chunk)))
        if index + MAX_LINE_LENGTH < len(nodes):
            open_file.write(",\n")
        else:
            open_file.write("\n")


class Geometry:  # pylint: disable=too-many-instance-attributes
    """Generates geometry for CalculiX simulations based on config file."""

    def __init__(self, config_path, output_path="solid.inp"):
        self.output_path = Path(output_path)
        self.cfg = read_config(config_path)

        self._all_nodes = []
        self._all_elements = []
        self._interface_faces = []
        self._current_node_id = 1
        self._current_element_id = 1
        self._fix_nodes = []
        self._wall_data = []
        self._interface_nodes = []
        self._interface_nodes_dict = None

        self.generate()

    def generate(self):
        """Generates geometry based on configuration."""
        orientations = [BOTTOM, TOP]
        self._interface_nodes_dict = {BOTTOM: [], TOP: []}
        for orientation in orientations:
            if self.cfg[WALLS][orientation]["active"]:
                self._wall_data.append(self.add_wall(orientation))

    def add_wall(self, orientation):  # noqa: WPS210
        """Adds a wall to the geometry.

        Args:
            orientation (str): Orientation of the wall ('bottom' or 'top').

        Returns:
            str: Wall String.
            int: Starting element ID of the wall.
            int: Ending element ID of the wall.
        """
        start_node_of_wall = self._current_node_id
        offset = self.cfg[WALLS][orientation]["offset"]
        nx = self.cfg[GEOMETRY]["n_elements_length"]
        nodes_per_row = nx + 1

        for row in range(2):
            self._generate_nodes(nodes_per_row, nx, offset, orientation, row)

        element_start = self._current_element_id
        for idx in range(nx):
            n1, n2, n3, n4 = self._add_vertices(
                idx, nodes_per_row, orientation, start_node_of_wall
            )

            self._all_elements.append(
                f"{self._current_element_id}, {n1}, {n2}, {n3}, {n4}"
            )
            self._current_element_id += 1

        face = "S3" if orientation == BOTTOM else "S1"
        self._interface_faces.append(f"E{orientation}, {face}")

        return f"E{orientation}", element_start, self._current_element_id - 1

    def write_file(  # noqa: WPS210, WPS213 pylint: disable=too-many-locals
            self, filepath, mesh_name="Nall", interface_name="Solid-Interface"
    ):
        """Writes the geometry to a CalculiX input file.

        Args:
            filepath (str or Path): Path to the output .inp file.
            mesh_name (str): Name of the mesh.
            interface_name (str): Name of the interface.

        Returns:
            Filepath without stem for simulation.
        """
        filepath = Path(filepath)
        work_dir = filepath.parent
        work_dir.mkdir(parents=True, exist_ok=True)

        inp_path = filepath.with_suffix(".inp")
        nam_filename = filepath.with_suffix(".nam").name
        nam_path = work_dir / nam_filename

        with open(nam_path, "w", encoding="utf-8") as nam_f:
            nam_f.write(f"*NSET, NSET=N{interface_name}\n")
            _write_node_list(nam_f, self._interface_nodes)

            nam_f.write(f"*SURFACE, NAME=S{interface_name}, TYPE=ELEMENT\n")
            for face_entry in self._interface_faces:
                nam_f.write(f"{face_entry}\n")

        material_config = self.cfg["material"]

        with open(inp_path, "w", encoding="uft-8") as inp_file:
            inp_file.write("** HEADING\n")
            inp_file.write(f"*NODE, NSET={mesh_name}\n")
            for node_str in self._all_nodes:
                inp_file.write(f"{node_str}\n")

            inp_file.write("*ELEMENT, TYPE=CPE4, ELSET=Eall\n")
            inp_file.write("\n".join(self._all_elements))
            inp_file.write("\n")

            for name, geo_set, element in self._wall_data:
                inp_file.write(
                    f"*ELSET, ELSET={name}, GENERATE\n{geo_set}, {element}, 1\n"
                )

            if self._fix_nodes:
                inp_file.write("*NSET, NSET=Nfix\n")
                _write_node_list(inp_file, self._fix_nodes)

            inp_file.write(f"*INCLUDE, INPUT={nam_filename}\n")

            inp_file.write("*MATERIAL, NAME=ELASTIC\n*ELASTIC\n")
            inp_file.write(
                f"{material_config["youngs_modulus"]}, "
                f"{material_config["poissons_ratio"]}\n"
            )
            inp_file.write(f"*DENSITY\n{material_config["density"]}\n")
            inp_file.write("*SOLID SECTION, ELSET=Eall, MATERIAL=ELASTIC\n1.0\n")

            inp_file.write("*STEP, NLGEOM, INC=1000000\n")
            inp_file.write("*DYNAMIC, DIRECT\n")
            dt = self.cfg.get("dt", DEFAULT_DT)
            duration = self.cfg.get("t_end", 1)
            inp_file.write(f"{dt}, {duration}\n")

            alpha = material_config.get("alpha", .0)
            beta = material_config.get("beta", DEFAULT_BETA)
            inp_file.write(f"*DAMPING, ALPHA={alpha}, BETA={beta}\n")

            inp_file.write("*BOUNDARY\n")
            if self._fix_nodes:
                inp_file.write("Nfix, 1, 2\n")
                inp_file.write("Nfix, 3, 3, 0.0\n")

            inp_file.write("*CLOAD\n")
            inp_file.write(f"N{interface_name}, 1, 0.0\n")
            inp_file.write(f"N{interface_name}, 2, 0.0\n")
            inp_file.write(f"N{interface_name}, 3, 0.0\n")

            inp_file.write(f"*NODE PRINT, NSET=N{interface_name}, FREQUENCY=1\n")
            inp_file.write("U, V, A\n")

            inp_file.write("*NODE FILE\nU\n")
            inp_file.write("*EL FILE\nS, E\n")
            inp_file.write("*END STEP\n")

        return filepath.parent / filepath.stem

    def _generate_nodes(  # pylint: disable=R0917, R0913
            self, nodes_per_row, nx, offset, orientation, row
    ):
        """Generates nodes for a wall.

        Args:
            nodes_per_row (int): Number of nodes per row.
            nx (int): Number of elements in length.
            offset (float): Offset of the wall.
            orientation (str): Orientation of the wall ('bottom' or 'top').
            row (int): Row index (0 or 1).
        """
        for index in range(nodes_per_row):
            self._add_row(index, nx, offset, orientation, row)

    def _add_row(    # pylint: disable=R0917, R0913
            self, index, nx, offset, orientation, row
    ):
        """Adds a row of nodes to the geometry.

        Args:
            index (int): Index of the node in the row.
            nx (int): Number of elements in length.
            offset (float): Offset of the wall.
            orientation (str): Orientation of the wall ('bottom' or 'top').
            row (int): Row index (0 or 1).
        """
        d_thick = row * self.cfg[GEOMETRY]["wall_thickness"]
        x_coord = index * (self.cfg[GEOMETRY]["length"] / nx)

        if orientation == BOTTOM:
            y_coord = offset - d_thick
        else:
            y_coord = (self.cfg[GEOMETRY]["height"] - offset) + d_thick

        self._all_nodes.append(
            f"{self._current_node_id}, {x_coord:.4f}, {y_coord:.4f}, 0.0"
        )

        if row == 0:
            self._interface_nodes_dict[orientation].append(self._current_node_id)
            self._interface_nodes.append(self._current_node_id)

        if self.cfg[WALLS][orientation].get("fixed_start") and index == 0:
            self._fix_nodes.append(self._current_node_id)
        if self.cfg[WALLS][orientation].get("fixed_end") and index == nx:
            self._fix_nodes.append(self._current_node_id)

        self._current_node_id += 1

    def _add_vertices(self, vertex_id, nodes_per_row, orientation, start_node_of_wall):
        """Calculates the vertex node IDs for an element.

        Args:
            vertex_id (int): Vertex index.
            nodes_per_row (int): Number of nodes per row.
            orientation (str): Orientation of the wall ('bottom' or 'top').
            start_node_of_wall (int): Starting node ID of the wall.

        Returns:
            tuple: Node IDs (n1, n2, n3, n4) of the element.
        """
        if orientation == BOTTOM:
            n1 = start_node_of_wall + nodes_per_row + vertex_id
            n2 = n1 + 1
            n3 = start_node_of_wall + vertex_id + 1
            n4 = n3 - 1
        else:
            n1 = start_node_of_wall + vertex_id
            n2 = n1 + 1
            n3 = n1 + nodes_per_row + 1
            n4 = n1 + nodes_per_row
        return n1, n2, n3, n4
