import pytest
from pymatgen.core import Structure, Lattice

from structure_builder.sublattice import get_sublattice_sites
from structure_builder.spec import CompositionSpec, SQSOptions
from structure_builder.sqs import generate_sqs


def tin_2x2x2():
    lat = Lattice.cubic(4.24)
    base = Structure(lat, ["Ti", "N"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    return base * [2, 2, 2]


def test_pymatgen_sqs_correct_composition():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = SQSOptions(iterations=200, backend="pymatgen")
    result = generate_sqs(s, sites, comp, opts)
    metal_species = [str(result[i].specie) for i in sites["metal"]]
    assert metal_species.count("Ti") == 4
    assert metal_species.count("Al") == 4


def test_auto_backend_works_without_atat():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = SQSOptions(iterations=200, backend="auto")
    result = generate_sqs(s, sites, comp, opts)
    assert len(result) == 16


def test_sqs_preserves_nonmetal_sublattice():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = SQSOptions(iterations=200, backend="pymatgen")
    result = generate_sqs(s, sites, comp, opts)
    nonmetal_species = [str(result[i].specie) for i in sites["nonmetal"]]
    assert all(sp == "N" for sp in nonmetal_species)


def test_sqs_nonmetal_substitution():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 1.0}, nonmetal={"N": 0.5, "C": 0.5})
    opts = SQSOptions(iterations=200, backend="pymatgen")
    result = generate_sqs(s, sites, comp, opts)
    nonmetal_species = [str(result[i].specie) for i in sites["nonmetal"]]
    assert nonmetal_species.count("N") == 4
    assert nonmetal_species.count("C") == 4


def test_unknown_backend_raises():
    s = tin_2x2x2()
    sites = get_sublattice_sites(s)
    comp = CompositionSpec(metal={"Ti": 0.5, "Al": 0.5})
    opts = SQSOptions(backend="bad_backend")
    with pytest.raises(ValueError, match="Unknown SQS backend"):
        generate_sqs(s, sites, comp, opts)
