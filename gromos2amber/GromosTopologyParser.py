from .gromos_format import parse_simple_columns, parse_array_block
from .Errors import GromosFormatError

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

    def getblock(self, blockname, checkheader = False):
        if not blockname in self.blocks:
            raise GromosFormatError(
                    "No block '"+blockname+"' in topology file."
            )
        block = self.blocks[blockname]
        if checkheader:
            try:
                numlines = int(block[1])
            except ValueError:
                message = "Could not parse the first line of block '{}' "\
                    " in the topology files as an integer: {}"
                raise GromosFormatError( message.format(blockname, block[1]))
            if not numlines+3 == len(block):
                message = "Expected {} lines in block '{}', found {} lines."
                raise GromosFormatError(
                    message.format(
                        numlines,
                        blockname,
                        len(block)-2,
                    )
                )

        return block

    def checkblockheader(self, blockname, size):
# for blocks where the first line is the number of data lines
        if blockname in self.blocks:
            raise GromosFormatError(
                    "No block '"+blockname+"' in topology file."
            )
        return self.blocks[blockname]

    def SOLUTEATOM(self):
        fieldwidths = [6,5,5,4,9,9,3,6]
        start_exclusions = sum(fieldwidths)
        start_neigh14_extra = sum(fieldwidths) - 2
        fieldbounds = [ ( sum(fieldwidths[0:i]), sum(fieldwidths[0:i+1]) ) 
                             for i in range(len(fieldwidths)) ]

        shortline = "The {}th line of the 'SOLUTEATOM' block in the topology "\
            "file is too short:\n'{}'"
        block = self.getblock("SOLUTEATOM")
        numatoms = int(block[1])
        atomindex, residue, typecode = [0]*numatoms, [0]*numatoms, [0]*numatoms
        name = [""]*numatoms
        mass, charge = [0.0]*numatoms, [0.0]*numatoms
        charge_group_code = [0]*numatoms
        exclusions = [ [] for a in range(numatoms) ]
        neigh14    = [ [] for a in range(numatoms) ]
        l = 2
        try:
            for i in range(numatoms):
                if len(block[l]) < sum(fieldwidths):
                    raise GromosFormatError(shortline.format(l, block[l]))

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
                    if len(block[l]) < start_exclusions + 1:
                        raise GromosFormatError(message.format(l, block[l]))
                    exclusions_string += block[l][start_exclusions:-1]
                exclusions[i] = [ int(exclusions_string[6*e:6*(e+1)])
                                  for e in range(num_exclusions) ]

                # assume space delimited (compatible with gromos 1.3.1 )
                l += 1
                neigh14_list =[ int(x) for x in block[l].split() ]
                num_neigh14 = neigh14_list[0] 
                neigh14[i].extend(neigh14_list[1:])
                while len(neigh14[i]) < num_neigh14:
                    l += 1
                    neigh14[i].extend([ int(x) for x in block[l].split() ])
                l += 1
        except ValueError as error:
            message = "Could not parse {}th atom on {}th "\
                "line of 'SOLUTEATOM' block. " + str(error)
            raise GromosFormatError(message.format(i, l))

        return ( atomindex, residue, name, typecode, mass, charge,
                 charge_group_code, exclusions, neigh14 )

    def TITLE(self):
        block = self.getblock("TITLE")
        return ''.join(block[1:])[0:-1] if len(block) > 1 else ''

    def PHYSICALCONSTANTS(self):
        block = self.getblock("PHYSICALCONSTANTS", checkheader = False)
        if 4+2 < len(block):
            raise GromosFormatError(
                "Only {} lines in PHYSICALCONSTANTS"\
                "block of topology file. Expected 4 lines."
            )
        return [ float(block[i+1]) for i in range(4) ]

    def TOPVERSION(self):
        return ''.join(self.getblock("TOPVERSION")[1])[0:-1]

    def ATOMTYPENAME(self):
        block = self.getblock("ATOMTYPENAME", checkheader = True)
        numtypes = int(block[1])
        return [ block[i+2][0:-1].strip() for i in range(numtypes) ]

    def RESNAMES(self):
        block = self.getblock("RESNAME", checkheader = True)
        numresidues = int(block[1])
        return [ block[i+2][0:-1] for i in range(numresidues) ]

    def BONDSTRETCHTYPE(self):
        block = self.getblock("BONDSTRETCHTYPE")
        return parse_simple_columns(block, [16,16,16], [float,float,float])

    def BOND(self, H = False):
        block = self.getblock( "BOND" + ("H" if H else "") )
        return parse_simple_columns(block, [7,7,5], [int,int,int])

    def BONDANGLEBENDTYPE(self):
        block = self.getblock("BONDANGLEBENDTYPE")
        return parse_simple_columns(block, [16,16,16], [float,float,float])

    def BONDANGLE(self, H = False):
        block = self.getblock( "BONDANGLE" + ("H" if H else "") )
        return parse_simple_columns(block, [7,7,7,5], [int,int,int,int])

    def IMPDIHEDRALTYPE(self):
        block = self.getblock("IMPDIHEDRALTYPE")
        return parse_simple_columns(block, [15,15], [float,float])

    def IMPDIHEDRAL(self, H = False):
        block = self.getblock( "IMPDIHEDRAL" + ("H" if H else "") )
        return parse_simple_columns(block, [7,7,7,7,5], [int,int,int,int,int])

    def TORSDIHEDRALTYPE(self):
        block = self.getblock("TORSDIHEDRALTYPE")
        return parse_simple_columns(block, [10,11,4], [float,float,int])

    def DIHEDRAL(self, H = False):
        block = self.getblock( "DIHEDRAL" + ("H" if H else "") )
        return parse_simple_columns(block, [7,7,7,7,5], [int,int,int,int,int])

    def LJPARAMETERS(self):
        block = self.getblock("LJPARAMETERS")
        return parse_simple_columns(block, [5,5,14,14,14,14],
                                    [int,int,float,float,float,float])
    def SOLUTEMOLECULES(self):
        block = self.getblock("SOLUTEMOLECULES")
        return parse_array_block(block, 6, int)

    def SOLVENTATOM(self):
        block = self.getblock("SOLVENTATOM")
        return parse_simple_columns(block, [4,6,4,11,11],
                                      [int,str,int,float,float])
    def SOLVENTCONSTR(self):
        block = self.getblock("SOLVENTCONSTR")
        return parse_simple_columns(block, [5,5,15], [int,int,float])

    def LJEXCEPTIONS(self):
        block = self.getblock("LJEXCEPTIONS")
        return parse_simple_columns(block, [5,5,14,14], [int,int,float,float] )


