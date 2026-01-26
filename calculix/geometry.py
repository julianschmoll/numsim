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
        nodes_per_row = 2 * nx + 1

        for row in range(3):
            d_thick = row * (self.cfg['geometry']['wall_thickness'] / 2.0)
            for index in range(nodes_per_row):
                x = index * (self.cfg["geometry"]["length"] / (2 * nx))

                if orientation == 'bottom':
                    y = offset - d_thick
                else:
                    y = (self.cfg["geometry"]["height"] - offset) + d_thick

                self._all_nodes.append(
                    f"{self._current_node_id}, {x:.8f}, {y:.8f}, 0.0" # Increased precision
                )

                if row == 0:
                    self._interface_nodes_dict[orientation].append(self._current_node_id)
                    self._interface_nodes.append(self._current_node_id) # Global list for Nall

                if self.cfg["walls"][orientation].get('fixed_start') and index == 0:
                    self._fix_nodes.append(self._current_node_id)
                if self.cfg["walls"][orientation].get('fixed_end') and index == 2 * nx:
                    self._fix_nodes.append(self._current_node_id)

                self._current_node_id += 1

        for i in range(nx):
            r0 = start_node_of_wall + (2 * i)
            r1 = r0 + nodes_per_row
            r2 = r1 + nodes_per_row

            if orientation == 'bottom':
                n1, n2, n3, n4 = r2, r2+2, r0+2, r0
                n5, n6, n7, n8 = r2+1, r1+2, r0+1, r1
            else:
                n1, n2, n3, n4 = r0, r0+2, r2+2, r2
                n5, n6, n7, n8 = r0+1, r1+2, r2+1, r1

            self._all_elements.append(
                f"{self._current_element_id}, "
                f"{n1}, {n2}, {n3}, {n4}, {n5}, {n6}, {n7}, {n8}"
            )
            self._current_element_id += 1

        face = "S3" if orientation == "bottom" else "S1"
        self._interface_faces.append(f"E{orientation}, {face}")

        return f"E{orientation}", self._current_element_id - nx, self._current_element_id - 1

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

            nam_f.write(f"*SURFACE, NAME=S{interface_name}\n")
            for face_entry in self._interface_faces:
                nam_f.write(f"{face_entry}\n")

        with open(inp_path, "w") as f:
            f.write("** HEADING\n")
            f.write(f"*NODE, NSET={mesh_name}\n")
            for node_str in self._all_nodes:
                f.write(f"{node_str}\n")

            f.write("*ELEMENT, TYPE=CPE8, ELSET=Eall\n")
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

            f.write("*STEP, INC=1000000\n")
            f.write("*DYNAMIC, DIRECT, NLGEOM\n")

            dt = self.cfg['geometry'].get('dt', 0.01)
            duration = self.cfg['geometry'].get('t_end', 10.0)
            f.write(f"{dt}, {duration}\n")
            f.write("*RESTART, WRITE, FREQUENCY=1\n")
            f.write("*BOUNDARY\n")
            if self._fix_nodes:
                f.write("Nfix, 1, 2\n")
            f.write(f"{mesh_name}, 3, 3\n")
            f.write("*CLOAD\n")
            f.write(f"N{interface_name}, 1, 0.0\n")
            f.write(f"N{interface_name}, 2, 0.0\n")

            f.write("*NODE FILE\nU\n")
            f.write("*EL FILE\nS, E\n")
            f.write("*END STEP\n")

            return filepath.parent / filepath.stem
