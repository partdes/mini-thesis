from __future__ import annotations
from typing import Dict, List

from pymatgen.core import Structure, Element


def get_sublattice_sites(structure: Structure) -> Dict[str, List[int]]:
    """
    Partition site indices into 'metal' and 'nonmetal' sublattices.
    Uses pymatgen's Element.is_metal and Element.is_metalloid for classification.
    """
    metal_sites: List[int] = []
    nonmetal_sites: List[int] = []

    for i, site in enumerate(structure):
        el = Element(str(site.specie).rstrip("0123456789+-"))
        if el.is_metal or el.is_metalloid:
            metal_sites.append(i)
        else:
            nonmetal_sites.append(i)

    return {"metal": metal_sites, "nonmetal": nonmetal_sites}
