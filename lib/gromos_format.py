
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
       raise Exception("No blocks parsed.")
   return blocks

def parse_simple_columns(block, widths, types, header = True):
    offset = 2 if header else 1
    nrows = len(block)-offset-1
    if header and nrows != int(block[1]):
        raise Exception("Wrong number of rows")
    ncols = len(widths)
    #check format is consistent with expectations
    if ncols != len(types):
        raise Exception("number of field widths not equal to number of types")
    line_width = sum(widths)+1 #includes newline
    for r in range(nrows):
        if len(block[r+offset]) != line_width:
            msg = "line {} of block {} is wrong length: \"{}\""
            raise Exception(msg.format(r+2, block[1].strip(), block[r+2]))
    # read columns into lists
    bounds = [ (sum(widths[0:i]) , sum(widths[0:i+1]))
                for i in range(len(widths)) ]
    return [
                [ types[c](block[i+offset][bounds[c][0]:bounds[c][1]])
                    for i in range(nrows) ]
            for c in range(ncols)
            ]

def parse_array_block(block, width, typ):
    n = int(block[1])
    line = ''.join(block[2:-1]).replace('\n','')
    if len(line) != n*width:
        msg = "Could not parse block:\n{}"
        raise Exception(msg.format(''.join(block)))
    return [ typ(line[i*width:(i+1)*width]) for i in range(n)]
