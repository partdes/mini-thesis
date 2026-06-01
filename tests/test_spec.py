import pytest
from pathlib import Path
from structure_builder.spec import load_spec, StructureSpec, CompositionSpec


def _write_spec(tmp_path, content):
    p = tmp_path / "spec.yaml"
    p.write_text(content, encoding="utf-8")
    return p


def test_load_minimal_spec(tmp_path):
    p = _write_spec(tmp_path, """
name: test_TiAlN
base_structure: data/structures/cubic_TiN.cif
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
""")
    spec = load_spec(p)
    assert spec.name == "test_TiAlN"
    assert spec.type == "compositional"
    assert spec.supercell_matrix == [2, 2, 2]
    assert spec.composition.metal == {"Ti": 0.5, "Al": 0.5}
    assert spec.composition.nonmetal == {}


def test_load_spec_with_nonmetal(tmp_path):
    p = _write_spec(tmp_path, """
name: test_TiAlCN
base_structure: data/structures/cubic_TiN.cif
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
  nonmetal:
    N: 0.75
    C: 0.25
""")
    spec = load_spec(p)
    assert spec.composition.nonmetal == {"N": 0.75, "C": 0.25}


def test_load_spec_sqs_options(tmp_path):
    p = _write_spec(tmp_path, """
name: test_sqs
base_structure: data/structures/cubic_TiN.cif
type: sqs
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
sqs_options:
  iterations: 1000
  backend: pymatgen
""")
    spec = load_spec(p)
    assert spec.sqs_options.iterations == 1000
    assert spec.sqs_options.backend == "pymatgen"


def test_load_spec_missing_name_raises(tmp_path):
    p = _write_spec(tmp_path, """
base_structure: data/structures/cubic_TiN.cif
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 1.0
""")
    with pytest.raises(ValueError, match="missing required fields"):
        load_spec(p)


def test_load_spec_invalid_type_raises(tmp_path):
    p = _write_spec(tmp_path, """
name: test
base_structure: data/structures/cubic_TiN.cif
type: bad_type
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 1.0
""")
    with pytest.raises(ValueError, match="type must be"):
        load_spec(p)


def test_load_spec_metal_fractions_not_1_raises(tmp_path):
    p = _write_spec(tmp_path, """
name: test
base_structure: data/structures/cubic_TiN.cif
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.3
    Al: 0.3
""")
    with pytest.raises(ValueError, match="sum to 1.0"):
        load_spec(p)


def test_load_spec_file_not_found():
    with pytest.raises(FileNotFoundError):
        load_spec("nonexistent/path/spec.yaml")


def test_load_spec_3x3_supercell_matrix(tmp_path):
    p = _write_spec(tmp_path, """
name: test_tilted
base_structure: data/structures/cubic_TiN.cif
type: compositional
supercell_matrix: [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
composition:
  metal:
    Ti: 1.0
""")
    spec = load_spec(p)
    assert spec.supercell_matrix == [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
