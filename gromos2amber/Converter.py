
from .Topology import Topology
from .Configuration import Configuration
from .AmberTopologyWriter import AmberTopologyWriter
from .AmberConfigurationWriter import AmberConfigurationWriter
from .Errors import GromosFormatError, IllegalArgumentError

def convert( topology_in,
                  topology_out,
                  config_in = None, 
                  config_out = None, 
                  solvent_resname="SOL",
                  num_solvent = -1,
                  ):
    if 4 < len(solvent_resname) and not 0 == len(solvent_resname):
        raise IllegalArgumentError(
            "Bad solvent residue name '{}'. ".format(solvent_resname) +\
                    "Solvent residue name must be 1-4 characters long."
        )
    
    if config_in == None and not config_out == None:
        raise IllegalArgumentError(
            "Output AMBER coordinates were requested but "\
                "no input gromos coordinates were provided."
        )
    
    try:
        topology = Topology(topology_in)
    except GromosFormatError as error:
        raise GromosFormatError( "Bad input topology format: " + str(error))
    
    if not config_in == None:
        try:
            config = Configuration(config_in)
            if not config.boxtype in (0, 1):
                raise GromosFormatError(
                    "Only vaccum or rectangular boxes are supported."
                )
        except GromosFormatError as error:
            raise GromosFormatError(
                "There is a problem with the coordinate file: " + str(error)
            )
        config.gather_molecules(topology)
        num_atoms = len(config.positions)
        num_solvent_molecules = ( num_atoms - len(topology.atoms) ) \
            * 1.0 / len(topology.solvent_atoms) 
        if not int(num_solvent_molecules) == num_solvent_molecules:
            raise GromosFormatError(
                "Mismatch between topology and coordinate files:"\
                "The apparent number of solvent atoms in the coordinate"\
                "file is {}. This is not divisible by the"\
                "number of atoms per solvent molecule ({})".format(
                    num_solute,
                    len(topology.solvent_atoms)
                )
            )

        num_solvent_molecules = int(num_solvent_molecules)
    else:
        num_solvent_molecules = num_solvent
    
    topology.add_solvent(num_solvent_molecules, solvent_resname)
    
    AmberTopologyWriter(topology).write(topology_out)
    
    if not config_out == None:
        AmberConfigurationWriter(config).write(config_out)

