from __future__ import annotations
import random
from typing import Any, Dict, List, Optional

from pymatgen.core import Structure, Element

from .spec import CompositionSpec, CompositionalOptions


def decorate(
    supercell: Structure,
    sublattice_sites: Dict[str, List[int]],
    composition: CompositionSpec,
    options: CompositionalOptions,
) -> Structure:
    """Apply compositional decoration to metal and optionally nonmetal sublattice."""
    result = _decorate_sublattice(
        supercell, sublattice_sites["metal"], composition.metal, options
    )
    if composition.nonmetal:
        nonmetal_opts = CompositionalOptions(arrangement=options.arrangement,
                                             layer_axis=options.layer_axis,
                                             site_occupancies=options.site_occupancies)
        result = _decorate_sublattice(
            result, sublattice_sites["nonmetal"], composition.nonmetal, nonmetal_opts
        )
    return result


def _decorate_sublattice(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
    options: CompositionalOptions,
) -> Structure:
    if options.arrangement == "random":
        return _random(structure, site_indices, composition)
    elif options.arrangement == "ordered":
        return _ordered(structure, site_indices, composition, options)
    elif options.arrangement == "explicit":
        return _explicit(structure, options.site_occupancies)
    else:
        raise ValueError(f"Unknown arrangement: {options.arrangement!r}")


def _random(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
) -> Structure:
    pool = _build_pool(site_indices, composition)
    random.shuffle(pool)
    new_species = [site.specie for site in structure]
    for idx, sp in zip(site_indices, pool):
        new_species[idx] = Element(sp)
    return Structure(structure.lattice, new_species, structure.frac_coords.copy())


def _ordered(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
    options: CompositionalOptions,
) -> Structure:
    ordered_type = options.ordered_type or "layered"
    if ordered_type in ("layered", "l10"):
        return _layered(structure, site_indices, composition, options.layer_axis)
    elif ordered_type == "l12":
        return _l12(structure, site_indices, composition)
    elif ordered_type == "b2":
        return _b2(structure, site_indices, composition)
    elif ordered_type == "checkerboard":
        return _checkerboard(structure, site_indices, composition)
    else:
        raise ValueError(f"Unknown ordered_type: {ordered_type!r}")


def _layered(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
    layer_axis: int,
) -> Structure:
    layers: Dict[float, List[int]] = {}
    for idx in site_indices:
        coord = round(structure[idx].frac_coords[layer_axis], 6)
        layers.setdefault(coord, []).append(idx)

    sorted_layers = sorted(layers.keys())
    elements = list(composition.keys())

    new_species = [site.specie for site in structure]
    for i, coord in enumerate(sorted_layers):
        el = elements[i % len(elements)]
        for idx in layers[coord]:
            new_species[idx] = Element(el)
    return Structure(structure.lattice, new_species, structure.frac_coords.copy())


def _l12(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
) -> Structure:
    """L12: minority species at face-center positions (two half-integer coords)."""
    elements = list(composition.keys())
    if len(elements) != 2:
        raise ValueError("l12 ordering requires exactly 2 species")
    majority, minority = elements[0], elements[1]

    new_species = [site.specie for site in structure]
    for idx in site_indices:
        fc = structure[idx].frac_coords
        half_int_count = sum(1 for c in fc if abs(c % 1 - 0.5) < 0.05)
        new_species[idx] = Element(minority if half_int_count == 2 else majority)
    return Structure(structure.lattice, new_species, structure.frac_coords.copy())


def _b2(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
) -> Structure:
    """B2: species A at corner-type sites, species B at body-center-type sites."""
    elements = list(composition.keys())
    if len(elements) != 2:
        raise ValueError("b2 ordering requires exactly 2 species")
    a_sp, b_sp = elements[0], elements[1]

    new_species = [site.specie for site in structure]
    for idx in site_indices:
        fc = structure[idx].frac_coords
        parity = round(sum(fc * 2)) % 2
        new_species[idx] = Element(a_sp if parity == 0 else b_sp)
    return Structure(structure.lattice, new_species, structure.frac_coords.copy())


def _checkerboard(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
) -> Structure:
    """Checkerboard: alternate in ab-plane by (h+k) parity."""
    elements = list(composition.keys())
    if len(elements) != 2:
        raise ValueError("checkerboard ordering requires exactly 2 species")
    a_sp, b_sp = elements[0], elements[1]

    new_species = [site.specie for site in structure]
    for idx in site_indices:
        fc = structure[idx].frac_coords
        parity = round(fc[0] * 4 + fc[1] * 4) % 2
        new_species[idx] = Element(a_sp if parity == 0 else b_sp)
    return Structure(structure.lattice, new_species, structure.frac_coords.copy())


def _explicit(
    structure: Structure,
    site_occupancies: Optional[List[Dict[str, Any]]],
) -> Structure:
    if not site_occupancies:
        raise ValueError("explicit arrangement requires a non-empty site_occupancies list")
    new_species = [site.specie for site in structure]
    for entry in site_occupancies:
        new_species[entry["index"]] = Element(entry["species"])
    return Structure(structure.lattice, new_species, structure.frac_coords.copy())


def _build_pool(site_indices: List[int], composition: Dict[str, float]) -> List[str]:
    n = len(site_indices)
    pool: List[str] = []
    remaining = n
    items = list(composition.items())
    for i, (element, fraction) in enumerate(items):
        count = remaining if i == len(items) - 1 else round(fraction * n)
        remaining -= count
        pool.extend([element] * count)
    return pool
