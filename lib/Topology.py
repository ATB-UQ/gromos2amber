
class Topology:

    def __init__(self, io):
        self.string = io.read()

    def write_amber(self, io):
        io.write(self.string)

