from fortran_format import fortran_format

class AmberConfigurationWriter:
    def __init__(self, configuration):
        self.configuration = configuration

    def write(self, io):
        title = "TEST TITLE"
        io.write(title+'\n')
        positions = self.configuration.positions
        io.write(fortran_format("I5,5E15.7",[len(positions),0]))
        values = []
        [ values.extend(pos) for pos in positions ]
        io.write(fortran_format("6F12.7", values))
        io.write(fortran_format("6F12.7", self.configuration.box_size))

