import pytest
from pymatgen.core import Structure, Lattice

from structure_builder.sublattice import get_sublattice_sites


def tin_primitive():
    lat = Lattice.cubic(4.24)
    return Structure(lat, ["Ti", "N"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def aln_primitive():
    lat = Lattice.cubic(4.05)
    return Structure(lat, ["Al", "N"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def test_tin_metal_site():
    s = tin_primitive()
    sites = get_sublattice_sites(s)
    species_at_metal = [str(s[i].specie) for i in sites["metal"]]
    assert "Ti" in species_at_metal
    assert "N" not in species_at_metal


def test_tin_nonmetal_site():
    s = tin_primitive()
    sites = get_sublattice_sites(s)
    species_at_nonmetal = [str(s[i].specie) for i in sites["nonmetal"]]
    assert "N" in species_at_nonmetal
    assert "Ti" not in species_at_nonmetal


def test_aln_partitions_correctly():
    s = aln_primitive()
    sites = get_sublattice_sites(s)
    assert len(sites["metal"]) == 1
    assert len(sites["nonmetal"]) == 1
    assert str(s[sites["metal"][0]].specie) == "Al"
    assert str(s[sites["nonmetal"][0]].specie) == "N"


def test_supercell_site_counts():
    s = tin_primitive() * [2, 2, 2]
    sites = get_sublattice_sites(s)
    assert len(sites["metal"]) == 8
    assert len(sites["nonmetal"]) == 8


def test_all_sites_covered():
    s = tin_primitive()
    sites = get_sublattice_sites(s)
    all_indices = sites["metal"] + sites["nonmetal"]
    assert sorted(all_indices) == list(range(len(s)))


def test_no_overlap_between_sublattices():
    s = tin_primitive() * [3, 3, 3]
    sites = get_sublattice_sites(s)
    assert set(sites["metal"]).isdisjoint(set(sites["nonmetal"]))
