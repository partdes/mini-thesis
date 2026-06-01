from __future__ import annotations
import os
from pathlib import Path
from typing import List

from pymatgen.core import Structure


def read_base_structure(source: str) -> Structure:
    """Read from local CIF/VASP path or fetch from Materials Project by mp-id."""
    if source.startswith("mp-"):
        return _fetch_from_mp(source)
    path = Path(source)
    if not path.exists():
        raise FileNotFoundError(f"Base structure not found: {path}")
    return Structure.from_file(str(path))


def write_structure(
    structure: Structure,
    name: str,
    formats: List[str],
    directory: str,
) -> List[Path]:
    """Write structure to each requested format. Returns list of written file paths."""
    out_dir = Path(directory)
    out_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    _SUPPORTED = {"vasp", "cif", "xyz"}

    for fmt in formats:
        fmt_lower = fmt.lower()
        if fmt_lower not in _SUPPORTED:
            raise ValueError(f"Unsupported output format: {fmt!r}. Choose from: vasp, cif, xyz")
        path = out_dir / f"{name}.{fmt_lower}"
        if fmt_lower == "vasp":
            structure.to(filename=str(path), fmt="poscar")
        elif fmt_lower == "cif":
            # Write string ourselves to avoid monty mode="w" vs "wt" incompatibility
            path.write_text(structure.to(fmt="cif"), encoding="utf-8")
        elif fmt_lower == "xyz":
            from pymatgen.io.xyz import XYZ
            XYZ(structure).write_file(str(path))
        written.append(path)

    return written


def _fetch_from_mp(mp_id: str) -> Structure:
    from dotenv import load_dotenv
    from mp_api.client import MPRester

    load_dotenv()
    load_dotenv("0.settings/.env")
    api_key = os.getenv("MP_API_KEY")
    if not api_key:
        raise EnvironmentError(
            "MP_API_KEY not set. Add MP_API_KEY=<key> to 0.settings/.env"
        )
    with MPRester(api_key) as mpr:
        return mpr.get_structure_by_material_id(mp_id)
