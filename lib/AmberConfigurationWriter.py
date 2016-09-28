from fortranformat import FortranRecordWriter as FortranWriter

DEBUG_VS_PDB = False

class AmberConfigurationWriter:
    def __init__(self, configuration):
        self.configuration = configuration

    def write(self, io):
        title = "TEST TITLE"
        io.write(title+'\n')
        positions = self.configuration.positions
        io.write(FortranWriter("I5,5E15.7").write([len(positions),0]) + '\n')
        values = []
        [ values.extend(pos) for pos in positions ]
        if DEBUG_VS_PDB:
            values = [ round(co, 3) for pos in values ]
        io.write(FortranWriter("6F12.7").write(values))
        io.write('\n')
        io.write(FortranWriter("6F12.7").write(self.configuration.box_size))
        io.write('\n')

