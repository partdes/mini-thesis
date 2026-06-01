import pytest
from collections import defaultdict
from pymatgen.core import Structure, Lattice

from structure_builder.sublattice import get_sublattice_sites
from structure_builder.spec import CompositionSpec, CompositionalOptions
from structure_builder.compositional import decorate


def tin_2x2x2():
    lat = Lattice.cubic(4.24)
    base = Structure(lat, ["Ti", "N"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    return base * [2, 2, 2]


def test_random_metal_species_counts():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = CompositionalOptions(arrangement="random")
    result = decorate(s, sites, comp, opts)
    metal_species = [str(result[i].specie) for i in sites["metal"]]
    assert metal_species.count("Ti") == 4
    assert metal_species.count("Al") == 4


def test_random_nonmetal_unchanged_when_single_species():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = CompositionalOptions(arrangement="random")
    result = decorate(s, sites, comp, opts)
    nonmetal_species = [str(result[i].specie) for i in sites["nonmetal"]]
    assert all(sp == "N" for sp in nonmetal_species)


def test_random_nonmetal_substituted_when_specified():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 1.0}, nonmetal={"N": 0.5, "C": 0.5})
    opts = CompositionalOptions(arrangement="random")
    result = decorate(s, sites, comp, opts)
    nonmetal_species = [str(result[i].specie) for i in sites["nonmetal"]]
    assert nonmetal_species.count("N") == 4
    assert nonmetal_species.count("C") == 4


def test_layered_each_layer_single_species():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = CompositionalOptions(arrangement="ordered", ordered_type="layered", layer_axis=2)
    result = decorate(s, sites, comp, opts)
    layers = defaultdict(set)
    for idx in sites["metal"]:
        layer = round(result[idx].frac_coords[2], 4)
        layers[layer].add(str(result[idx].specie))
    for layer_species in layers.values():
        assert len(layer_species) == 1, f"Layer has mixed species: {layer_species}"


def test_l12_correct_site_assignment():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.75, "Al": 0.25})
    opts = CompositionalOptions(arrangement="ordered", ordered_type="l12")
    result = decorate(s, sites, comp, opts)
    metal_species = [str(result[i].specie) for i in sites["metal"]]
    assert "Ti" in metal_species
    assert "Al" in metal_species


def test_explicit_assigns_specified_sites():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    metal_sites = sites["metal"]
    occupancies = [
        {"index": metal_sites[0], "species": "Al"},
        {"index": metal_sites[1], "species": "Ti"},
    ]
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = CompositionalOptions(arrangement="explicit", site_occupancies=occupancies)
    result = decorate(s, sites, comp, opts)
    assert str(result[metal_sites[0]].specie) == "Al"
    assert str(result[metal_sites[1]].specie) == "Ti"


def test_explicit_missing_occupancies_raises():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = CompositionalOptions(arrangement="explicit", site_occupancies=None)
    with pytest.raises(ValueError, match="site_occupancies"):
        decorate(s, sites, comp, opts)


def test_unknown_arrangement_raises():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 1.0})
    opts = CompositionalOptions(arrangement="unknown_type")
    with pytest.raises(ValueError, match="Unknown arrangement"):
        decorate(s, sites, comp, opts)
