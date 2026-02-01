from ccx2paraview import Converter
from pathlib import Path
import shutil


def convert_to_vtk(file_to_convert, output_path):
    file_to_convert = Path(file_to_convert)

    output_path = Path(output_path)
    cleaned_suffix = output_path.suffix.split(".")[-1]

    if cleaned_suffix not in ["vtk", "vtu"]:
        raise RuntimeError("Unsupported file type requested")

    converter = Converter(str(file_to_convert), [cleaned_suffix])
    converter.run()

    # The converter does not let us specify an output filename or path
    output_path = Path(output_path)
    converter_output = (file_to_convert.parent
                       / f"{file_to_convert.stem}.{cleaned_suffix}")

    if converter_output == output_path:
        return output_path

    if not output_path.parent.exists():
        output_path.parent.mkdir(parents=True, exist_ok=True)

    vtk_files = sorted(file_to_convert.parent.glob(f"{file_to_convert.stem}*.{cleaned_suffix}"))

    for idx, vtk_file in enumerate(vtk_files):
        shutil.copy(
            vtk_file,
            output_path.parent / f"{output_path.stem}_{idx:04d}.{cleaned_suffix}"
        )
    return output_path
