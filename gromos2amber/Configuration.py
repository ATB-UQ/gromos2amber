from . import gromos_format as gf
import sys

NANOMETRE = 10.0
PICOSECONDS = 20.455 #amber time unit (1/20.455 ps)

class Configuration:
    def __init__(self, io):
        blocks = gf.parse_blocks(io)
        nm = NANOMETRE
        ps = PICOSECONDS
        if "GENBOX" in blocks:
            genbox = blocks["GENBOX"]
            self.box_size = [ float(x)*nm for x in genbox[2].split() ]
        elif "BOX" in blocks:
            box = blocks["BOX"]
            self.box_size = [ float(x)*nm for x in box[1].split() ]
        else:
            self.box_size = [0.0, ] * 3

        posblock, cols, types  = (
            "POSITION",
            [5, 6, 6, 7, 15, 15, 15],
            (int, str, str, int, float, float, float),
        ) if "POSITION" in blocks else (
            "POSITIONRED",
            [15, 15, 15],
            (float, float, float),
        ) if "POSITIONRED" in blocks else (None, None, None)

        if None == posblock:
            raise GromosFormatError(
                "No 'POSITION' or 'POSITIONRED' block found "\
                    "in coordinate file"
            )
        columns = gf.parse_simple_columns(
            blocks[posblock],
            cols,
            types,
            header = False,
        )

        x,y,z = columns[-3:]

        if "LATTICESHIFTS" in blocks:
            sx,sy,sz = gf.parse_simple_columns(blocks["LATTICESHIFTS"],
                                               [10,10,10],
                                               (int,int,int),
                                               header = False)
        else:
            numatoms = len(x)
            sx = [0]*numatoms
            sy,sz = list(sx), list(sx)

        if "VELOCITY" in blocks or "VELOCITYRED" in blocks:
            velblock, cols, types = ("VELOCITY",
                                     [5,6,6,7,15,15,15],
                                     (int,str,str,int,float,float,float),
                                     ) if "VELOCITY" in blocks \
                                             else ("VELOCITYRED",
                                                   [15,15,15],
                                                   (float,float,float),
                                                   )
            columns = gf.parse_simple_columns(blocks[velblock],
                                              cols,
                                              types,
                                              header = False)
            vx,vy,vz = columns[-3:]

            self.velocities = [ [xi*nm/ps, yi*nm/ps, zi*nm/ps]
                                    for xi,yi,zi in zip(vx,vy,vz) ]
        else:
            self.velocities = None

        bx,by,bz = self.box_size

        self.positions = [ [xi*nm+sxi*bx, yi*nm+syi*by, zi*nm+szi*bz]
                            for xi,yi,zi,sxi,syi,szi in zip(x,y,z,sx,sy,sz) ]
        self.title = ''.join(blocks["TITLE"][1:-1]).strip()

    def gather_molecules(self, topology):
        x = self.positions
        box = self.box_size
        bond_lists = (topology.bonds_wH, topology.bonds_woH)
        if sum(box) == 0:
            return []
        num_broken_bond_dims = 1
        while num_broken_bond_dims > 0:
            num_broken_bond_dims = 0
            for bonds in bond_lists:
                for bond in bonds:
                    i,j = bond.atoms
                    for d in range(3):
                        if abs(x[i][d]-x[j][d]) > 0.5*box[d]:
                            num_broken_bond_dims += 1
                            if x[i][d] > x[j][d]:
                                if x[j][d]>box[d]: raise(Exception("stuck in loop"))
                                x[j][d] += box[d]
                            else:
                                if x[i][d]>box[d]: raise(Exception("stuck in loop"))
                                x[i][d] += box[d]

