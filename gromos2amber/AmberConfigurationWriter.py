from .fortran_format import fortran_format

class AmberConfigurationWriter:
    def __init__(self, configuration):
        self.configuration = configuration

    def write(self, io):
        title = self.configuration.title.replace('\n','; ')
        io.write(title+'\n')
        positions = self.configuration.positions
        io.write(fortran_format("i5,5e15.7",[len(positions),0]))
        xvalues = []
        [ xvalues.extend(pos) for pos in positions ]
        io.write(fortran_format("6f12.7", xvalues))
        velocities = self.configuration.velocities
        if not velocities == None:
            vvalues = []
            [ vvalues.extend(vel) for vel in velocities ]
            io.write(fortran_format("6f12.7", vvalues))
        io.write(fortran_format("6f12.7", self.configuration.box_size))

