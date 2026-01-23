import os
from pathlib import Path


def get_template(template_name):
    """Helper to get template content."""
    template_path = (
            Path(__file__).parent.parent / "resources" / "openfoam"
            / f"{template_name}.txt"
    )
    with open(template_path, 'r') as f:
        return f.read()


def write_file(path, content):
    """Helper to write content to a file."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        f.write(content)


def generate_header(class_name, object_name, location=None):
    header_template = get_template("header")
    loc_str = f'    location    "{location}";\n' if location else ""
    return header_template.format(
        class_name=class_name, object_name=object_name, loc_str=loc_str
    )


def boundary_condition(bc_name, value_x, value_y, is_pressure=False, config=None, side=None):
    """Maps custom BC names to OpenFOAM types."""
    is_pulsatile = (side == 'left' and not is_pressure and config.get('amplitude', 0) != 0)

    if bc_name == "InflowNoSlip":
        if is_pressure:
            return "type zeroGradient;"
        if is_pulsatile:
            return get_template("pulsatile_inlet").format(
                u_mean_x=value_x, u_mean_y=value_y, amp=config['amplitude'],
                freq=config['frequency'], shift=config['timeShift'],
            )
        return f"type fixedValue; value uniform ({value_x} {value_y} 0);"

    if bc_name == "Outflow" and is_pressure:
        return "type fixedValue; value uniform 0;"

    return "type zeroGradient;"


def generate_u(c):
    u_template = get_template("U")
    return generate_header("volVectorField", "U") + u_template.format(
        left_bc=boundary_condition(
            c['boundaryLeft'], c['dirichletLeftX'], c['dirichletLeftY'],
            False, c, 'left'
        ),
        right_bc=boundary_condition(
            c['boundaryRight'], c['dirichletRightX'], c['dirichletRightY'],
            False, c, 'right'
        ),
        bottom_bc=boundary_condition(
            c['boundaryBottom'], c['dirichletBottomX'], c['dirichletBottomY'],
            False, c,'bottom'
        ),
        top_bc=boundary_condition(
            c['boundaryTop'],c['dirichletTopX'], c['dirichletTopY'],
            False, c, 'top'),
    )


def generate_p(c):
    p_template = get_template("p")
    return generate_header(
        "volScalarField", "p"
    ) + p_template.format(
        left_bc=boundary_condition(
            c['boundaryLeft'], c['dirichletLeftX'], c['dirichletLeftY'],
            True, c, 'left'
        ),
        right_bc=boundary_condition(
            c['boundaryRight'], c['dirichletRightX'], c['dirichletRightY'],
            True, c, 'right'
        ),
        bottom_bc=boundary_condition(
            c['boundaryBottom'], c['dirichletBottomX'], c['dirichletBottomY'],
            True, c,'bottom'
        ),
        top_bc=boundary_condition(
            c['boundaryTop'], c['dirichletTopX'], c['dirichletTopY'],
            True, c, 'top'
        )
    )

def generate_g(c):
    return ""

def generate_point_displacement(c):
    return generate_header(
        "pointVectorField", "pointDisplacement"
    ) + get_template("pointDisplacement")


def generate_block_mesh(c):
    z_thick = 0.1
    blockmesh_template = get_template("blockMesh")
    return generate_header(
        "dictionary", "blockMeshDict"
    ) + blockmesh_template.format(
        Lx=c['physicalSizeX'],Ly=c['physicalSizeY'], Lz=z_thick,
        nCellsX=c['nCellsX'], nCellsY=c['nCellsY'], nCellsZ=1
    )


def generate_control_dict(c):
    tpl = get_template("controlDict")
    return generate_header(
        "dictionary", "controlDict"
    ) + tpl.format(endTime=c['endTime'], maximumDt=c['maximumDt'])


def generate_fv_schemes(c):
    div_scheme = "Gauss linearUpwind grad(U)" if c.get('useDonorCell') else "Gauss linear"
    tpl = get_template("fvSchemes")
    return generate_header(
        "dictionary", "fvSchemes"
    ) + tpl.format(div_scheme=div_scheme)


def generate_fv_solution(c):
    tpl = get_template("fvSolution")
    return generate_header(
        "dictionary", "fvSolution"
    ) + tpl.format()


def generate_decompose_par_dict(c):
    tpl = get_template("decomposeParDict")
    return generate_header(
        "dictionary", "decomposeParDict"
    ) + tpl.format(nSubdomains=c.get('nSubdomains', 4))


def generate_precice_dict(c):
    interface_patch = c.get('fsiInterface', "top")
    tpl = get_template("preciceDict")
    return generate_header(
        "dictionary", "preciceDict"
    ) + tpl.format(interface_patch=interface_patch)


def generate_transport_properties(c):
    L = float(c['physicalSizeY'])
    U_char = c.get('dirichletLeftX', 1.0)
    if U_char == 0:
        U_char = 1.0
    Re = float(c['re'])
    nu = (U_char * L) / Re
    tpl = get_template("transportProperties")
    return generate_header(
        "dictionary", "transportProperties"
    ) + tpl.format(nu=f"{nu:.6e}")


def generate_turbulence_properties(c):
    tpl = get_template("turbulenceProperties")
    return generate_header(
        "dictionary", "turbulenceProperties"
    ) + tpl.format()


def generate_dynamic_mesh_dict(c):
    tpl = get_template("dynamicMeshDict")
    return generate_header(
        "dictionary", "dynamicMeshDict"
    ) + tpl.format()


def get_bc(bc_name, val_x, val_y, field_type, config, side):
    """
    field_type: 'U', 'p', or 'pointDisplacement'

    `config['fsiInterface']` may be:
      - a string: "top" or "top,bottom"
      - a list/tuple/set: ["top", "bottom"]
    """
    raw = config.get('fsiInterface', "top,bottom")
    if isinstance(raw, (list, tuple, set)):
        fsi_interfaces = {str(x).strip().lower() for x in raw}
    else:
        fsi_interfaces = {p.strip().lower() for p in str(raw).split(",") if p.strip()}

    side_lc = str(side).strip().lower()

    if field_type == 'pointDisplacement':
        if side_lc == 'frontandback' or side_lc == 'frontandback':
            return "type empty;"
        if side_lc in fsi_interfaces:
            return "type fixedValue; value uniform (0 0 0);"
        return "type fixedValue; value uniform (0 0 0);"

    if field_type == 'p':
        if side_lc == 'frontandback':
            return "type empty;"
        if bc_name == "Outflow":
            return "type fixedValue; value uniform 0;"
        return "type zeroGradient;"

    if field_type == 'U':
        if side_lc == 'frontandback':
            return "type empty;"
        if side_lc in fsi_interfaces:
            return "type movingWallVelocity; value uniform (0 0 0);"

        if bc_name == "InflowNoSlip":
            if side_lc == 'left' and config.get('amplitude', 0) != 0:
                amp, freq, shift = config['amplitude'], config['frequency'], config['timeShift']
                return f"""
type            codedFixedValue;
value           uniform ({val_x} {val_y} 0);
name            pulsatileInlet;
code
#{{
    const scalar t = this->db().time().value();
    const scalar amp = {amp};
    const scalar freq = {freq};
    const scalar shift = {shift};
    const scalar pi = constant::mathematical::pi;
    scalar scale = 1.0 + amp * sin(2*pi*freq*t + shift);
    operator== (vector({val_x}, {val_y}, 0) * scale);
#}};
                """
            else:
                return f"type fixedValue; value uniform ({val_x} {val_y} 0);"

        elif bc_name == "Outflow":
            return "type zeroGradient;"

    return "type zeroGradient;"



def generate_field(c, name):
    cls = "volScalarField" if name == 'p' else "volVectorField"
    if name == 'pointDisplacement':
        cls = "pointVectorField"

    dims = "[0 2 -2 0 0 0 0]" if name == 'p' else "[0 1 0 0 0 0 0]" if name == 'pointDisplacement' else "[0 1 -1 0 0 0 0]"

    default_internal = "0" if name == 'p' else "(0 0 0)"

    content = generate_header(cls, name) + f"""
dimensions      {dims};
internalField   uniform {default_internal};

boundaryField
{{
"""
    for side in ['left', 'right', 'bottom', 'top', 'frontAndBack']:
        bc_key = f"boundary{side.capitalize()}"
        val_x_key = f"dirichlet{side.capitalize()}X"
        val_y_key = f"dirichlet{side.capitalize()}Y"

        bc_name = c.get(bc_key, "NoSlip")  # Default if missing
        val_x = c.get(val_x_key, 0)
        val_y = c.get(val_y_key, 0)

        content += f"    {side}\n    {{\n        {get_bc(bc_name, val_x, val_y, name, c, side)}\n    }}\n"

    content += "}\n"

    return content
