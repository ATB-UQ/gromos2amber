
def parse_blocks(io):
   blocks = {}
   lines = [line for line in io.readlines()
               if line.strip() and not line[0] == '#' ]
   start = 0
   for l,line in enumerate(lines):
       if l == start:
           blockname = line.strip()
       elif line == "END\n":
           blocks[blockname] = lines[start:l+1]
           start = l+1
       else:
           continue
   if len(blocks) == 0:
       raise FormatException("No blocks found in file.")
   return blocks

def parse_simple_columns(block, widths, types, header = True):
    offset = 2 if header else 1
    nrows = len(block)-offset-1
    if header:
        header_nrows = int(block[1])
        if nrows != header_nrows:
            message = "Block '{}' contains {} lines of data, expected {}"
            raise FormatException(
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
            raise FormatException(
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
        raise FormatError(
            "Block '{}' could not be parsed".format(block[0].strip())
        )
    return columns

def parse_array_block(block, width, typ):
    n = int(block[1])
    line = ''.join(block[2:-1]).replace('\n','')
    if len(line) != n*width:
        raise FormatError(
            "Block '{}' could not be parsed".format(block[0].strip())
        )
    try:
        result = [ typ(line[i*width:(i+1)*width]) for i in range(n)]
    except ValueError as error:
        raise FormatError(
            "Block '{}' could not be parsed. ".format(block[0].strip())+\
                str(error)
        )
    return result
