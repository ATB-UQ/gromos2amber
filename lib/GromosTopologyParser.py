
class GromosTopologyParser:

    def __init__(self, io):
        self.blocks = {}
        lines = [line for line in io.readlines()
                    if line.strip() and not line[0] == '#' ]
        start = 0
        for l,line in enumerate(lines):
            if l == start:
                blockname = line.strip()
            elif line == "END\n":
                self.blocks[blockname] = lines[start:l+1]
                start = l+1
            else:
                continue
        if len(self.blocks) == 0:
            raise Exception("No blocks parsed.")

    def SOLUTEATOM(self):
        fieldwidths = [6,5,5,4,9,9,3,6]
        start_exclusions = sum(fieldwidths)
        fieldbounds = [ ( sum(fieldwidths[0:i]), sum(fieldwidths[0:i+1]) ) 
                             for i in range(len(fieldwidths)) ]
        block = self.blocks["SOLUTEATOM"]
        numatoms = int(block[1])
        atomindex, residue, typecode = [0]*numatoms, [0]*numatoms, [0]*numatoms
        name = [""]*numatoms
        mass, charge = [0.0]*numatoms, [0.0]*numatoms
        charge_group_code = [0]*numatoms
        exclusions = [ [] for a in range(numatoms) ]
        neigh14    = [ [] for a in range(numatoms) ]
        l = 2
        for i in range(numatoms):
            fields = [block[l][a:b] for a,b in fieldbounds ]
            atomindex[i] = int(fields[0])
            residue[i]   = int(fields[1])
            name[i]      = fields[2].strip()
            typecode[i]  = int(fields[3])
            mass[i], charge[i] = float(fields[4]), float(fields[5])
            charge_group_code[i] = int(fields[6])
            num_exclusions = int(fields[7])

            exclusions_string = block[l][sum(fieldwidths):-1]
            while len(exclusions_string) < 6*num_exclusions:
                l += 1
                exclusions_string += block[l][start_exclusions:-1]
            exclusions[i] = [ int(exclusions_string[6*e:6*(e+1)])
                              for e in range(num_exclusions) ]

            l += 1
            num_neigh14 = int(block[l][0:start_exclusions])
            neigh14_string = block[l][sum(fieldwidths):-1]
            while len(neigh14_string) < 6*num_neigh14:
                l += 1
                neigh14_string += block[l][start_exclusions:-1]
            neigh14[i] = [ int(neigh14_string[6*n:6*(n+1)])
                              for n in range(num_neigh14) ]
            l += 1

        return ( atomindex, residue, name, typecode, mass, charge,
                 charge_group_code, exclusions, neigh14 )

    def TITLE(self):
        return ''.join(self.blocks["TITLE"][1:])[0:-1]

    def PHYSICALCONSTANTS(self):
        block = self.blocks["PHYSICALCONSTANTS"]
        return [ float(block[i+1]) for i in range(4) ]

    def TOPVERSION(self):
        return ''.join(self.blocks["TOPVERSION"][1])[0:-1]

    def ATOMTYPENAME(self):
        block = self.blocks["ATOMTYPENAME"]
        numtypes = int(block[1])
        return [ block[i+2][0:-1].strip() for i in range(numtypes) ]

    def RESNAMES(self):
        block = self.blocks["RESNAME"]
        numresidues = int(block[1])
        return [ block[i+2][0:-1] for i in range(numresidues) ]

    def BONDSTRETCHTYPE(self):
        block = self.blocks["BONDSTRETCHTYPE"]
        return parse_simple_columns(block, [16,16,16], [float,float,float])

    def BOND(self, H = False):
        block = self.blocks[ "BOND" + ("H" if H else "") ]
        return parse_simple_columns(block, [7,7,5], [int,int,int])

    def BONDANGLEBENDTYPE(self):
        block = self.blocks["BONDANGLEBENDTYPE"]
        return parse_simple_columns(block, [16,16,16], [float,float,float])

    def BONDANGLE(self, H = False):
        block = self.blocks[ "BONDANGLE" + ("H" if H else "") ]
        return parse_simple_columns(block, [7,7,7,5], [int,int,int,int])

    def IMPDIHEDRALTYPE(self):
        block = self.blocks["IMPDIHEDRALTYPE"]
        return parse_simple_columns(block, [15,15], [float,float])

    def IMPDIHEDRAL(self, H = False):
        block = self.blocks[ "IMPDIHEDRAL" + ("H" if H else "") ]
        return parse_simple_columns(block, [7,7,7,7,5], [int,int,int,int,int])

    def TORSDIHEDRALTYPE(self):
        block = self.blocks["TORSDIHEDRALTYPE"]
        return parse_simple_columns(block, [10,11,4], [float,float,int])

    def DIHEDRAL(self, H = False):
        block = self.blocks[ "DIHEDRAL" + ("H" if H else "") ]
        return parse_simple_columns(block, [7,7,7,7,5], [int,int,int,int,int])

    def LJPARAMETERS(self):
        block = self.blocks["LJPARAMETERS"]
        return parse_simple_columns(block, [5,5,14,14,14,14],
                                    [int,int,float,float,float,float])
    def SOLUTEMOLECULES(self):
        block = self.blocks["SOLUTEMOLECULES"]
        parse_array_block(block, 6, int)

    def SOLVENTATOM(self):
        block = self.blocks["SOLVENTATOM"]
        return parse_simple_columns(block, [4,6,4,11,11],
                                      [int,str,int,float,float])

    def LJEXCEPTIONS(self):
        block = self.blocks["LJEXCEPTIONS"]
        return parse_simple_columns(block, [5,5,14,14], [int,int,float,float] )


    
def parse_simple_columns(block, widths, types):
    nrows = int(block[1])
    ncols = len(widths)
    #check format is consistent with expectations
    if ncols != len(types):
        raise Exception("number of field widths not equal to number of types")
    line_width = sum(widths)+1 #includes newline
    for r in range(nrows):
        if len(block[r+2]) != line_width:
            msg = "line {} of block {} is wrong length: \"{}\""
            raise Exception(msg.format(r+2, block[1].strip(), block[r+2]))
    # read columns into lists
    bounds = [ (sum(widths[0:i]) , sum(widths[0:i+1]))
                for i in range(len(widths)) ]
    return [
      [ types[c](block[i+2][bounds[c][0]:bounds[c][1]]) for i in range(nrows) ]
      for c in range(ncols)
    ]

def parse_array_block(block, width, typ):
    n = int(block[1])
    line = ''.join(block[2:]).replace('\n','')
    if len(line) != n*width:
        raise Exception("Could not parse {}".format(block[0]))
    return [ typ(line[i*n:(i+1)*n]) for i in range(n)]



