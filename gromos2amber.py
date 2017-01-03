
import sys
from lib.Topology import Topology
from lib.Configuration import Configuration
from lib.AmberTopologyWriter import AmberTopologyWriter
from lib.AmberConfigurationWriter import AmberConfigurationWriter

def gromos2amber( topology_in,
                  topology_out,
                  config_in = None, 
                  config_out = None, 
                  solvent_resname="SOL",
                  num_solvent = -1,
                  ):
    if len(solvent_resname)>4:
        raise(Exception(
                "ERROR: Solvent residue name cannot be longer than 4 characers."))
    
    if solvent_resname == "":
        raise(Exception(
                "ERROR: Solvent residue name cannot be an empty string."))
    
    if config_in == None and not config_out == None:
        sys.stderr.write("WARNING: Cannot write configuration file when no input "
                         +"configuration file has been supplied.")
        config_out = None
    
    topology = Topology(topology_in)
    
    if not config_in == None:
        config = Configuration(config_in)
        config.gather_molecules(topology)
        num_atoms = len(config.positions)
        num_solvent_molecules = ( num_atoms - len(topology.atoms) ) \
                                * 1.0 / len(topology.solvent_atoms) 
        if not int(num_solvent_molecules) == num_solvent_molecules:
            msg = "Number of solvent atoms ({}) not divisible by number of "\
                    "number of atoms per solvent molecule ({})"
            raise Exception(msg.format(num_solute,
                                       len(topology.solvent_atoms)))
        num_solvent_molecules = int(num_solvent_molecules)
    else:
        num_solvent_molecules = num_solvent
    
    topology.add_solvent(num_solvent_molecules, solvent_resname)
    
    AmberTopologyWriter(topology).write(topology_out)
    
    if not config_out == None:
        AmberConfigurationWriter(config).write(config_out)
    

