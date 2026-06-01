import pytest
from pathlib import Path

from structure_builder.builder import build_from_spec

TIN_CIF = Path("data/structures/cubic_TiN.cif")


def _write_spec(tmp_path, content):
    p = tmp_path / "spec.yaml"
    p.write_text(content, encoding="utf-8")
    return str(p)


@pytest.mark.skipif(not TIN_CIF.exists(), reason="requires cubic_TiN.cif")
def test_build_compositional_random(tmp_path):
    spec_path = _write_spec(tmp_path, f"""
name: test_TiAlN
base_structure: {TIN_CIF}
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
output:
  formats: [vasp]
  directory: {tmp_path}
compositional_options:
  arrangement: random
""")
    written = build_from_spec(spec_path)
    assert (tmp_path / "test_TiAlN.vasp").exists()
    assert len(written) == 1


@pytest.mark.skipif(not TIN_CIF.exists(), reason="requires cubic_TiN.cif")
def test_build_sqs(tmp_path):
    spec_path = _write_spec(tmp_path, f"""
name: test_TiAlN_sqs
base_structure: {TIN_CIF}
type: sqs
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
output:
  formats: [vasp, cif]
  directory: {tmp_path}
sqs_options:
  iterations: 100
  backend: pymatgen
""")
    written = build_from_spec(spec_path)
    assert (tmp_path / "test_TiAlN_sqs.vasp").exists()
    assert (tmp_path / "test_TiAlN_sqs.cif").exists()
    assert len(written) == 2


@pytest.mark.skipif(not TIN_CIF.exists(), reason="requires cubic_TiN.cif")
def test_build_layered(tmp_path):
    spec_path = _write_spec(tmp_path, f"""
name: test_TiAlN_layered
base_structure: {TIN_CIF}
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 0.5
    Al: 0.5
output:
  formats: [vasp]
  directory: {tmp_path}
compositional_options:
  arrangement: ordered
  ordered_type: layered
  layer_axis: 2
""")
    written = build_from_spec(spec_path)
    assert (tmp_path / "test_TiAlN_layered.vasp").exists()


def test_build_unknown_type_raises(tmp_path):
    spec_path = _write_spec(tmp_path, f"""
name: bad
base_structure: {TIN_CIF}
type: compositional
supercell_matrix: [2, 2, 2]
composition:
  metal:
    Ti: 1.0
output:
  formats: [vasp]
  directory: {tmp_path}
compositional_options:
  arrangement: random
""")
    # This should succeed (Ti only, no substitution needed)
    written = build_from_spec(spec_path)
    assert len(written) == 1
