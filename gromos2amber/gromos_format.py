from .Errors import GromosFormatError
from chemistry_helpers.gromos_blocks import GromosBlockError, blocks_dict


def parse_blocks(io, empty_message="No blocks found in file."):
   """{name: lines}, each list framed as the callers index it: name line, value lines,
   "END\\n" (every line newline-terminated, fixed-width lines kept unstripped)."""
   try:
       split = blocks_dict(io.read(), on_duplicate="last", strip=False)
   except GromosBlockError as error:
       raise GromosFormatError(str(error))
   if not split:
       raise GromosFormatError(empty_message)
   return {name: [name + "\n"] + [line + "\n" for line in lines] + ["END\n"]
           for name, lines in split.items()}

def parse_simple_columns(block, widths, types, header = True):
    offset = 2 if header else 1
    nrows = len(block)-offset-1
    if header:
        header_nrows = int(block[1])
        if nrows != header_nrows:
            message = "Block '{}' contains {} lines of data, expected {}"
            raise GromosFormatError(
                message.format(block[0].strip()),
                nrows,
                header_nrows,
            )
    ncols = len(widths)
    #check format is consistent with expectations
    if ncols != len(types):
        raise Exception(
            "Number of field widths not equal to number of types"
        )
    line_width = sum(widths)+1 #includes newline
    for r in range(nrows):
        if len(block[r+offset]) != line_width:
            message = "line {} of block is wrong length: \"{}\""
            raise GromosFormatError(
                message.format(r+2, block[1].strip(), block[r+2])
            )
    # read columns into lists
    bounds = [ (sum(widths[0:i]) , sum(widths[0:i+1]))
                for i in range(len(widths)) ]
    try:
        columns = [
            [ types[c](block[i+offset][bounds[c][0]:bounds[c][1]])
                for i in range(nrows) ]
            for c in range(ncols)
        ]
    except ValueError:
        raise GromosFormatError(
            "Block '{}' could not be parsed".format(block[0].strip())
        )
    return columns

def parse_array_block(block, width, typ):
    n = int(block[1])
    line = ''.join(block[2:-1]).replace('\n','')
    if len(line) != n*width:
        raise GromosFormatError(
            "Block '{}' could not be parsed".format(block[0].strip())
        )
    try:
        result = [ typ(line[i*width:(i+1)*width]) for i in range(n)]
    except ValueError as error:
        raise GromosFormatError(
            "Block '{}' could not be parsed. ".format(block[0].strip())+\
                str(error)
        )
    return result
