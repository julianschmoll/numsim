# python
from pathlib import Path
import subprocess
import shutil

def export_vtk(case_dir, out_dir=None):
    """
    Run foamToVTK and copy the latest time VTK/VTU/PVTu files into out_dir.
    Returns the output directory Path.
    """
    case_dir = Path(case_dir).resolve()
    subprocess.run(["foamToVTK", "-case", str(case_dir)], check=True, cwd=case_dir)

    vtk_root = case_dir / "VTK"
    if not vtk_root.exists():
        # fallback: search for any .vtk/.vtu files under case_dir
        search_root = case_dir
    else:
        # choose latest numeric time subdirectory under VTK if present
        times = [d for d in vtk_root.iterdir() if d.is_dir() and all(ch.isdigit() or ch == '.' for ch in d.name)]
        if times:
            search_root = max(times, key=lambda p: float(p.name))
        else:
            search_root = vtk_root

    out_dir = Path(out_dir or (case_dir.parent / "out")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    for ext in ("*.vtk", "*.vtu", "*.pvtu"):
        for f in search_root.rglob(ext):
            shutil.copy2(f, out_dir / f.name)

    return out_dir
