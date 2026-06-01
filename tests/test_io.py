import pytest
from pathlib import Path
from pymatgen.core import Structure, Lattice

from structure_builder.io import read_base_structure, write_structure

TIN_CIF = Path("data/structures/cubic_TiN.cif")


def simple_tin():
    lat = Lattice.cubic(4.24)
    return Structure(lat, ["Ti", "N"], [[0, 0, 0], [0.5, 0.5, 0.5]])


def test_read_local_cif():
    if not TIN_CIF.exists():
        pytest.skip("cubic_TiN.cif not found")
    structure = read_base_structure(str(TIN_CIF))
    species = {str(s.specie) for s in structure}
    assert "Ti" in species
    assert "N" in species


def test_read_local_vasp():
    vasp_path = Path("data/structures/cubic_TiN.vasp")
    if not vasp_path.exists():
        pytest.skip("cubic_TiN.vasp not found")
    structure = read_base_structure(str(vasp_path))
    species = {str(s.specie) for s in structure}
    assert "Ti" in species


def test_read_nonexistent_raises():
    with pytest.raises(FileNotFoundError, match="not found"):
        read_base_structure("data/structures/does_not_exist.cif")


def test_write_vasp(tmp_path):
    s = simple_tin()
    written = write_structure(s, "test_tin", ["vasp"], str(tmp_path))
    assert (tmp_path / "test_tin.vasp").exists()
    assert len(written) == 1


def test_write_cif(tmp_path):
    s = simple_tin()
    written = write_structure(s, "test_tin", ["cif"], str(tmp_path))
    assert (tmp_path / "test_tin.cif").exists()


def test_write_xyz(tmp_path):
    s = simple_tin()
    written = write_structure(s, "test_tin", ["xyz"], str(tmp_path))
    assert (tmp_path / "test_tin.xyz").exists()


def test_write_multiple_formats(tmp_path):
    s = simple_tin()
    written = write_structure(s, "test_tin", ["vasp", "cif"], str(tmp_path))
    assert len(written) == 2
    assert (tmp_path / "test_tin.vasp").exists()
    assert (tmp_path / "test_tin.cif").exists()


def test_write_creates_directory(tmp_path):
    s = simple_tin()
    out_dir = tmp_path / "new_subdir"
    write_structure(s, "test_tin", ["vasp"], str(out_dir))
    assert (out_dir / "test_tin.vasp").exists()


def test_write_unsupported_format_raises(tmp_path):
    s = simple_tin()
    with pytest.raises(ValueError, match="Unsupported output format"):
        write_structure(s, "test", ["pdb"], str(tmp_path))
