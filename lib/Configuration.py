from . import gromos_format as gf

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

        _,_,_,_,x,y,z = gf.parse_simple_columns(blocks["POSITION"],
                                            [5,6,6,7,15,15,15],
                                            (int,str,str,int,float,float,float),
                                            header = False)

        if "LATTICESHIFTS" in blocks:
            sx,sy,sz = gf.parse_simple_columns(blocks["LATTICESHIFTS"],
                                               [10,10,10],
                                               (int,int,int),
                                               header = False)
        else:
            numatoms = len(x)
            sx = [0]*numatoms
            sy,sz = list(sx), list(sx)

        if "VELOCITY" in blocks:
            _,_,_,_,vx,vy,vz = gf.parse_simple_columns(blocks["VELOCITY"],
                                            [5,6,6,7,15,15,15],
                                            (int,str,str,int,float,float,float),
                                            header = False)
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
                                x[j][d] += box[d]
                            else:
                                x[i][d] += box[d]

