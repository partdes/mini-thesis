from __future__ import annotations
from pathlib import Path
from typing import List

from pymatgen.core import Structure

from .spec import load_spec
from .io import read_base_structure, write_structure
from .sublattice import get_sublattice_sites
from .compositional import decorate
from .sqs import generate_sqs


def build_from_spec(spec_path: str | Path) -> List[Path]:
    """
    Build a supercell from a YAML spec file.
    Returns list of written output file paths.
    """
    spec = load_spec(spec_path)
    base = read_base_structure(spec.base_structure)
    supercell = base * spec.supercell_matrix
    sublattice_sites = get_sublattice_sites(supercell)

    if spec.type == "compositional":
        from .spec import CompositionalOptions
        opts = spec.compositional_options or CompositionalOptions()
        structure = decorate(supercell, sublattice_sites, spec.composition, opts)
    elif spec.type == "sqs":
        from .spec import SQSOptions
        opts = spec.sqs_options or SQSOptions()
        structure = generate_sqs(supercell, sublattice_sites, spec.composition, opts)
    else:
        raise ValueError(f"Unknown type: {spec.type!r}")

    return write_structure(
        structure,
        spec.name,
        spec.output.formats,
        spec.output.directory,
    )
