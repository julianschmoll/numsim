import yaml
from pathlib import Path


def read_config(config_path):
    if not Path(config_path).exists():
        raise FileNotFoundError
    if not Path(config_path).is_file():
        raise FileNotFoundError
    with open(config_path, "r") as config_file:
        return yaml.load(config_file, Loader=yaml.FullLoader)


class Geometry:
    def __init__(self, config_path, output_path="solid.inp"):
        self.output_path = Path(output_path)
        try:
            self.cfg = read_config(config_path)
        except Exception as error:
            raise RuntimeError("Could not read config file") from error

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
        orientations = ["bottom", "top"]
        self._interface_nodes_dict = {"bottom": [], "top": []}
        for orientation in orientations:
            if self.cfg['walls'][orientation]['active']:
                self._wall_data.append(self.add_wall(orientation))

    def add_wall(self, orientation):
        start_node_of_wall = self._current_node_id
        offset = self.cfg["walls"][orientation]["offset"]
        nx = self.cfg['geometry']['n_elements_length']
        nodes_per_row = nx + 1

        # Generate nodes
        for row in range(2):
            d_thick = row * self.cfg['geometry']['wall_thickness']
            for index in range(nodes_per_row):
                x = index * (self.cfg["geometry"]["length"] / nx)

                if orientation == 'bottom':
                    # row 0 is interface (top), row 1 is base (bottom)
                    y = offset - d_thick
                else:
                    # row 0 is interface (bottom), row 1 is base (top)
                    y = (self.cfg["geometry"]["height"] - offset) + d_thick

                self._all_nodes.append(f"{self._current_node_id}, {x:.4f}, {y:.4f}, 0.0")

                if row == 0:
                    self._interface_nodes_dict[orientation].append(self._current_node_id)
                    self._interface_nodes.append(self._current_node_id)

                if self.cfg["walls"][orientation].get('fixed_start') and index == 0:
                    self._fix_nodes.append(self._current_node_id)
                if self.cfg["walls"][orientation].get('fixed_end') and index == nx:
                    self._fix_nodes.append(self._current_node_id)

                self._current_node_id += 1

        # Generate quadrilateral elements with CCW logic
        element_start = self._current_element_id
        for i in range(nx):
            if orientation == "bottom":
                # For bottom wall, Row 1 is physically BELOW Row 0
                n1 = start_node_of_wall + nodes_per_row + i  # R1 bottom-left
                n2 = n1 + 1                                  # R1 bottom-right
                n3 = start_node_of_wall + i + 1              # R0 top-right
                n4 = n3 - 1                                  # R0 top-left
            else:
                # For top wall, Row 0 is physically BELOW Row 1
                n1 = start_node_of_wall + i                  # R0 bottom-left
                n2 = n1 + 1                                  # R0 bottom-right
                n3 = n1 + nodes_per_row + 1                  # R1 top-right
                n4 = n1 + nodes_per_row                      # R1 top-left

            self._all_elements.append(f"{self._current_element_id}, {n1}, {n2}, {n3}, {n4}")
            self._current_element_id += 1

        face = "S3" if orientation == "bottom" else "S1"
        self._interface_faces.append(f"E{orientation}, {face}")

        return f"E{orientation}", element_start, self._current_element_id - 1
    def _write_node_list(self, fileobj, nodes):
        """Helper to write a list of nodes formatted for CalculiX"""
        for i in range(0, len(nodes), 16):
            chunk = nodes[i:i+16]
            fileobj.write(", ".join(map(str, chunk)))
            if i + 16 < len(nodes):
                fileobj.write(",\n")
            else:
                fileobj.write("\n")

    def write_file(self, filepath, mesh_name="Nall", interface_name="Solid-Interface"):
        filepath = Path(filepath)
        work_dir = filepath.parent
        work_dir.mkdir(parents=True, exist_ok=True)

        inp_path = filepath.with_suffix(".inp")
        nam_filename = filepath.with_suffix(".nam").name
        nam_path = work_dir / nam_filename

        with open(nam_path, "w") as nam_f:
            nam_f.write(f"*NSET, NSET=N{interface_name}\n")
            self._write_node_list(nam_f, self._interface_nodes)

            nam_f.write(f"*SURFACE, NAME=S{interface_name}, TYPE=ELEMENT\n")
            for face_entry in self._interface_faces:
                nam_f.write(f"{face_entry}\n")

        with open(inp_path, "w") as f:
            f.write("** HEADING\n")
            f.write(f"*NODE, NSET={mesh_name}\n")
            for node_str in self._all_nodes:
                f.write(f"{node_str}\n")

            f.write("*ELEMENT, TYPE=CPE4, ELSET=Eall\n")
            f.write("\n".join(self._all_elements) + "\n")

            for name, s, e in self._wall_data:
                f.write(f"*ELSET, ELSET={name}, GENERATE\n{s}, {e}, 1\n")

            if self._fix_nodes:
                f.write("*NSET, NSET=Nfix\n")
                self._write_node_list(f, self._fix_nodes)

            f.write(f"*INCLUDE, INPUT={nam_filename}\n")

            f.write("*MATERIAL, NAME=ELASTIC\n*ELASTIC\n")
            f.write(f"{self.cfg['material']['youngs_modulus']}, {self.cfg['material']['poissons_ratio']}\n")
            f.write(f"*DENSITY\n{self.cfg['material']['density']}\n")
            f.write("*SOLID SECTION, ELSET=Eall, MATERIAL=ELASTIC\n1.0\n")

            f.write("*STEP, NLGEOM, INC=1000000\n")
            f.write("*DYNAMIC\n")

            alpha = 0.0  # Mass proportional damping
            beta = 0.001  # Stiffness proportional damping
            f.write(f"*DAMPING, ALPHA={alpha}, BETA={beta}\n")

            dt = self.cfg['geometry'].get('dt', 0.01)
            duration = self.cfg['geometry'].get('t_end', 10.0)
            f.write(f"{dt}, {duration}\n")

            f.write("*BOUNDARY\n")
            if self._fix_nodes:
                f.write("Nfix, 1, 2\n")
                f.write("Nfix, 3, 3, 0.0\n")

            f.write("*CLOAD\n")
            f.write(f"N{interface_name}, 1, 0.0\n")
            f.write(f"N{interface_name}, 2, 0.0\n")
            f.write(f"N{interface_name}, 3, 0.0\n")

            f.write(f"*NODE PRINT, NSET=N{interface_name}, FREQUENCY=1\n")
            f.write("U, V, A\n")  # Displacement, velocity, acceleration


            f.write("*NODE FILE\nU\n")
            f.write("*EL FILE\nS, E\n")
            f.write("*END STEP\n")

            return filepath.parent / filepath.stem
