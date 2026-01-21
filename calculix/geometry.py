import yaml

def generate(config_path="2d_elastic_tube.yaml"):
    try:
        with open(config_path, 'r') as f:
            cfg = yaml.safe_load(f)
    except FileNotFoundError:
        print("Error: solid_settings.yaml not found.")
        return

    # Extract Physics/Geometry
    L = cfg['geometry']['length_x']
    H = cfg['geometry']['length_y']
    thick = cfg['geometry']['wall_thickness']
    nx = cfg['geometry']['n_elements_length']

    all_nodes = []
    all_elements = []
    fix_nodes = []
    interface_faces = []

    current_node_id = 1
    current_elem_id = 1

    def add_wall_logic(orientation):
        nonlocal current_node_id, current_elem_id
        start_node_of_wall = current_node_id
        wall_cfg = cfg['walls'][orientation]
        offset = wall_cfg['offset']
        nodes_per_row = 2 * nx + 1

        # 1. Generate Nodes
        for row in range(3):
            d_thick = row * (thick / 2.0)
            for i in range(nodes_per_row):
                x = i * (L / (2 * nx))
                if orientation == 'bottom':
                    y = offset - d_thick
                else:
                    y = (H - offset) + d_thick

                all_nodes.append(f"{current_node_id}, {x:.4f}, {y:.4f}, 0.0")

                # Granular Boundary Logic
                if wall_cfg.get('fixed_start') and i == 0:
                    fix_nodes.append(current_node_id)
                if wall_cfg.get('fixed_end') and i == 2 * nx:
                    fix_nodes.append(current_node_id)

                current_node_id += 1

        # 2. Elements (CPE8) - Strictly Counter-Clockwise
        for i in range(nx):
            r0 = start_node_of_wall + (2 * i)      # Inner
            r1 = r0 + nodes_per_row                # Middle
            r2 = r1 + nodes_per_row                # Outer

            if orientation == 'bottom':
                n1, n2, n3, n4 = r2, r2+2, r0+2, r0
                n5, n6, n7, n8 = r2+1, r1+2, r0+1, r1
            else:
                n1, n2, n3, n4 = r0, r0+2, r2+2, r2
                n5, n6, n7, n8 = r0+1, r1+2, r2+1, r1

            all_elements.append(f"{current_elem_id}, {n1}, {n2}, {n3}, {n4}, {n5}, {n6}, {n7}, {n8}")
            current_elem_id += 1

        face = "S3" if orientation == "bottom" else "S1"
        interface_faces.append(f"E{orientation}, {face}")
        return f"E{orientation}", current_elem_id - nx, current_elem_id - 1

    # Execute
    elsets = []
    for orient in ['bottom', 'top']:
        if cfg['walls'][orient]['active']:
            elsets.append(add_wall_logic(orient))

    # Write solid.inp
    with open("solid.inp", "w") as f:
        f.write("** HEADING\n*NODE, NSET=Nall\n")
        f.write("\n".join(all_nodes) + "\n")
        f.write("*ELEMENT, TYPE=CPE8, ELSET=Eall\n")
        f.write("\n".join(all_elements) + "\n")
        for name, s, e in elsets:
            f.write(f"*ELSET, ELSET={name}, GENERATE\n{s}, {e}, 1\n")
        if fix_nodes:
            f.write("*NSET, NSET=Nfix\n" + ", ".join(map(str, fix_nodes)) + "\n")
        f.write("*MATERIAL, NAME=ELASTIC\n*ELASTIC\n")
        f.write(f"{cfg['material']['youngs_modulus']}, {cfg['material']['poissons_ratio']}\n")
        f.write(f"*DENSITY\n{cfg['material']['density']}\n")
        f.write("*SOLID SECTION, ELSET=Eall, MATERIAL=ELASTIC\n1.0\n")
        f.write("*SURFACE, NAME=Sinterface\n" + "\n".join(interface_faces) + "\n")
        f.write("*STEP, NLGEOM\n*STATIC\n*BOUNDARY\n")
        if fix_nodes: f.write("Nfix, 1, 3\n")
        f.write("Nall, 3, 3\n")
        f.write("*NODE FILE\nU\n*EL FILE\nS\n*END STEP\n")
    print("Generated solid.inp")

    return "solid.inp"


if __name__ == "__main__":
    generate()
