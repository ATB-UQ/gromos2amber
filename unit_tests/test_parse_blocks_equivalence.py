"""The shared splitter behind parse_blocks agrees with the pre-migration copy on well-formed input."""
import io
import pytest
from gromos2amber.gromos_format import parse_blocks
from gromos2amber.Errors import GromosFormatError


def old_parse_blocks(f):  # verbatim pre-migration copy
    blocks = {}
    lines = [line for line in f.readlines() if line.strip() and not line[0] == '#']
    start = 0
    for l, line in enumerate(lines):
        if l == start:
            blockname = line.strip()
        elif line == "END\n":
            blocks[blockname] = lines[start:l + 1]
            start = l + 1
    if len(blocks) == 0:
        raise GromosFormatError("No blocks found in file.")
    return blocks


TOP = ("# comment\nTITLE\nmol # not a comment\nEND\n\nSOLUTEATOM\n    2\n"
       "     1    1 C     12.0110    0.0000  0   1    2\n     2    1 H      1.0080    0.0000  0   0\nEND\n"
       "BOX\n  1.0 2.0 3.0\nEND\n")
CASES = [TOP, "TITLE\nx\nEND\nBOX\n1 2 3\nEND\n", "BOX\n1 2 3\nEND\nBOX\n4 5 6\nEND\n"]


@pytest.mark.parametrize("text", CASES)
def test_same_as_old_copy(text):
    assert parse_blocks(io.StringIO(text)) == old_parse_blocks(io.StringIO(text))


def test_empty_message_is_kept():
    with pytest.raises(GromosFormatError, match="custom"):
        parse_blocks(io.StringIO(""), "custom")


def test_unterminated_block_is_now_an_error():
    with pytest.raises(GromosFormatError):
        parse_blocks(io.StringIO("BOX\n1 2 3\n"))
