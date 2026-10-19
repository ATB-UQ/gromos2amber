"""Unit tests for gromos2amber. Run with the worktree on PYTHONPATH."""
import io
import os

import pytest

from gromos2amber import GromosFormatError, IllegalArgumentError, convert
from gromos2amber.fortran_format import fortran_format
from gromos2amber.gromos_format import parse_blocks, parse_simple_columns
from gromos2amber.Topology import (
    Atom, Interaction, KILOJOULE, NANOMETRE, DEGREE,
    _canonicalise_exclusions, _create_solvent,
)

BENCH_TOP = "/home/atb/ATB/ncmas_benchmarks/systems/ti_8040/bench.top"


def test_unit_constants():
    assert KILOJOULE == pytest.approx(1 / 4.184)
    assert NANOMETRE == 10.0
    assert DEGREE == pytest.approx(0.017453292519943295)


def test_fortran_format_wraps_and_handles_empty():
    assert fortran_format("10i8", list(range(12))) == (
        "".join("{:>8d}".format(i) for i in range(10)) + "\n"
        + "{:>8d}{:>8d}".format(10, 11) + "\n")
    assert fortran_format("10i8", []) == "\n"


def test_parse_blocks_splits_on_END():
    blocks = parse_blocks(io.StringIO("# c\nTITLE\nhello\nEND\nBOX\n1 2 3\nEND\n"))
    assert set(blocks) == {"TITLE", "BOX"}


def test_parse_blocks_empty_file_raises():
    with pytest.raises(GromosFormatError):
        parse_blocks(io.StringIO(""))


def test_parse_simple_columns_reads_fixed_widths():
    block = ["X\n", "2\n", "  1  2\n", "  3  4\n", "END\n"]
    assert parse_simple_columns(block, [3, 3], [int, int]) == [[1, 3], [2, 4]]


@pytest.mark.xfail(strict=True, raises=IndexError,
                   reason="review finding: error message built with too few format args")
def test_parse_simple_columns_header_count_mismatch_is_a_format_error():
    block = ["X\n", "3\n", "  1  2\n", "  3  4\n", "END\n"]
    with pytest.raises(GromosFormatError):
        parse_simple_columns(block, [3, 3], [int, int])


def test_exclusions_stored_once_on_lower_atom():
    atoms = [Atom("A", 0, 1.0, 0.0, [1, 2], []),
             Atom("B", 0, 1.0, 0.0, [0], []),
             Atom("C", 0, 1.0, 0.0, [], [])]
    atoms[2].exclusions = [0]          # recorded only on the higher-indexed partner
    _canonicalise_exclusions(atoms)
    assert [a.exclusions for a in atoms] == [[1, 2], [], []]


def test_solvent_exclusions_point_forward_only():
    sol = [Atom("OW", 0, 16.0, -0.8, [], []), Atom("H1", 1, 1.0, 0.4, [], []),
           Atom("H2", 1, 1.0, 0.4, [], [])]
    bonds = [Interaction([0, 1], 0), Interaction([0, 2], 0)]
    atoms, b = _create_solvent(sol, bonds, 2, 10)
    assert [a.exclusions for a in atoms[:3]] == [[11, 12], [12], []]
    assert [a.exclusions for a in atoms[3:]] == [[14, 15], [15], []]
    assert b[2].atoms == [13, 14]


def test_bad_solvent_resname_rejected():
    with pytest.raises(IllegalArgumentError):
        convert(io.StringIO(""), io.StringIO(), solvent_resname="TOOLONG")


def test_coordinates_out_without_in_rejected():
    with pytest.raises(IllegalArgumentError):
        convert(io.StringIO(""), io.StringIO(), config_out=io.StringIO())


@pytest.mark.skipif(not os.path.exists(BENCH_TOP), reason="ncmas bench.top not present")
def test_real_topology_converts_with_expected_pointers():
    out = io.StringIO()
    with open(BENCH_TOP) as fh:
        convert(fh, out, num_solvent=0)
    text = out.getvalue()
    assert text.startswith("%VERSION")
    lines = text.splitlines()
    i = lines.index("%FLAG POINTERS")
    assert int(lines[i + 2].split()[0]) == 44      # natom of the 44-atom solute
    # GROMOS CHB (kJ mol-1 nm-2) -> kcal mol-1 A-2 and the AMBER 1/2 convention:
    k_first_bond = float(lines[lines.index("%FLAG BOND_FORCE_CONSTANT") + 2].split()[0])
    assert 100 < k_first_bond < 1500
