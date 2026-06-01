from __future__ import annotations
import random
from typing import Dict, List, Tuple

from pymatgen.core import Structure, Element

from .spec import CompositionSpec, SQSOptions
from .compositional import _random


def generate_sqs(
    supercell: Structure,
    sublattice_sites: Dict[str, List[int]],
    composition: CompositionSpec,
    options: SQSOptions,
) -> Structure:
    """Generate an SQS structure using the requested backend."""
    if options.backend == "atat":
        return _atat_sqs(supercell, sublattice_sites, composition, options)
    elif options.backend == "pymatgen":
        return _mc_sqs(supercell, sublattice_sites, composition, options.iterations)
    elif options.backend == "auto":
        try:
            return _atat_sqs(supercell, sublattice_sites, composition, options)
        except (NotImplementedError, Exception):
            return _mc_sqs(supercell, sublattice_sites, composition, options.iterations)
    else:
        raise ValueError(f"Unknown SQS backend: {options.backend!r}. Choose: atat, pymatgen, auto")


def _mc_sqs(
    supercell: Structure,
    sublattice_sites: Dict[str, List[int]],
    composition: CompositionSpec,
    iterations: int,
) -> Structure:
    """Monte Carlo SQS optimizer minimizing Warren-Cowley first-shell SRO parameters."""
    metal_sites = sublattice_sites["metal"]
    nonmetal_sites = sublattice_sites["nonmetal"]

    current = _random(supercell, metal_sites, composition.metal)
    if composition.nonmetal:
        current = _random(current, nonmetal_sites, composition.nonmetal)

    best = current.copy()
    best_score = _wc_score(current, metal_sites, composition.metal)

    for _ in range(iterations):
        candidate, swapped = _try_swap(current, metal_sites)
        if not swapped:
            continue
        score = _wc_score(candidate, metal_sites, composition.metal)
        if score < best_score:
            best = candidate.copy()
            best_score = score
            current = candidate

    return best


def _wc_score(
    structure: Structure,
    site_indices: List[int],
    composition: Dict[str, float],
) -> float:
    """
    Sum of squared Warren-Cowley parameters for the first coordination shell.
    Lower score = more random-like pair correlations.
    """
    cutoff = 3.5
    site_set = set(site_indices)
    species_at = {i: str(structure[i].specie) for i in site_indices}

    pair_counts: Dict[Tuple[str, str], int] = {}
    bond_counts: Dict[str, int] = {}

    for i in site_indices:
        sp_i = species_at[i]
        for neighbor in structure.get_neighbors(structure[i], cutoff):
            j = neighbor.index
            if j not in site_set or j == i:
                continue
            sp_j = species_at[j]
            pair_counts[(sp_i, sp_j)] = pair_counts.get((sp_i, sp_j), 0) + 1
            bond_counts[sp_i] = bond_counts.get(sp_i, 0) + 1

    score = 0.0
    for (sp_i, sp_j), count in pair_counts.items():
        total = bond_counts.get(sp_i, 1)
        actual = count / total
        ideal = composition.get(sp_j, 0.0)
        score += (actual - ideal) ** 2
    return score


def _try_swap(
    structure: Structure, site_indices: List[int]
) -> Tuple[Structure, bool]:
    """Swap two randomly chosen sites with different species. Returns (structure, did_swap)."""
    for _ in range(50):
        i, j = random.sample(site_indices, 2)
        sp_i = str(structure[i].specie)
        sp_j = str(structure[j].specie)
        if sp_i != sp_j:
            new_species = [site.specie for site in structure]
            new_species[i] = Element(sp_j)
            new_species[j] = Element(sp_i)
            return Structure(structure.lattice, new_species, structure.frac_coords.copy()), True
    return structure, False


def _atat_sqs(supercell, sublattice_sites, composition, options):
    raise NotImplementedError(
        "ATAT backend requires ATAT installed separately. "
        "Use backend='pymatgen' or backend='auto' instead."
    )
