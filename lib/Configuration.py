import gromos_format as gf

class Configuration:

    def __init__(self, io):
        blocks = gf.parse_blocks(io)
        block = blocks["POSITION"]
        _,_,_,_,x,y,z = gf.parse_simple_columns(block, [5,6,6,7,15,15,15],
                                            (int,str,str,int,float,float,float),
                                            header = False)
        self.positions = [ [xi,yi,zi] for xi,yi,zi in zip(x,y,z) ]


