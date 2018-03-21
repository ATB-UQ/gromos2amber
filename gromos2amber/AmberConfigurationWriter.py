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
        if 0 != sum(self.configuration.box_angle):
            box = self.configuration.box_size + self.configuration.box_angle
        else:
            box = self.configuration.box_size

        io.write(fortran_format("6f12.7", box))


