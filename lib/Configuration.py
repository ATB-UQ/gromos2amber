import gromos_format as gf

NANOMETRE = 10.0

class Configuration:
    def __init__(self, io):
        blocks = gf.parse_blocks(io)
        nm = NANOMETRE
        genbox = blocks["GENBOX"]
        self.box_size = [ float(x)*nm for x in genbox[2].split() ]

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

        bx,by,bz = self.box_size

        self.positions = [ [xi*nm+sxi*bx, yi*nm+syi*by, zi*nm+szi*bz]
                            for xi,yi,zi,sxi,syi,szi in zip(x,y,z,sx,sy,sz) ]



