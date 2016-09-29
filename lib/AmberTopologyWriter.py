""" Functions which return information required to generate a parm7 file.

Each function is named after the corresponding section, taking ATB data as
the only argument, and returning 4 values:

  1. A list of values for the body of the section
  2. The Fortran format string
  3. A section comment to appear after the "%FLAG ..." line
  4. An integer used to sort the sections relative to each other

if no comment is required for the section (because it is self-explantitory,
eg "TITLE"), the amber_sections.NOCOMMENT is returned as the comment.

  """

# Because amber.py uses the inspect module to assign functions in
# amber_sections.py to sections in the self.topology file, it is important
# that the functions in amber_sections all correspond to sections. Helper
# functions have been put in amber_helpers.py. Unfortunately, we cannot simply
# import all the amber_helpers functions, as inspect will then treat them as
# though they are part of amber_sections.py.

from fortranformat import FortranRecordWriter as FortranWriter
from inspect import getmembers, ismethod

# Constants
NOCOMMENT = "" # comment when no comment required
NOTERMS = "This section intentionally left empty."
LINEWIDTH = 80

class AmberTopologyWriter:

    def __init__(self, topology):
        self.topology = topology

    def write(self, io):
        not_sections = ["write", "__init__"]
        section_functions = [ member
                            for member in getmembers(self,
                                                     predicate=ismethod)
                            if not member[0] in not_sections]
        sections = []
        ordering = []

        for title,func in section_functions:
            values, format_string, comment, order = func()
            sections.append( _section(title, comment, format_string, values) )
            ordering.append( order )

        io.write("%VERSION  VERSION_STAMP = V0001.000\n")

        for order, sec in sorted(zip(ordering, sections)):
            io.write(sec)

    def CTITLE(self):
        format_string = '20A4'
        comment = NOCOMMENT
        order = 0
        title = self.topology.get_title()
        title = title + " "*(LINEWIDTH - len(title))
        values = [title[i:i+4] for i in range(0, LINEWIDTH, 4)]
        return values, format_string, comment, order
    
    def FORCE_FIELD_TYPE(self):
        format_string ='I2,A78'
        comment = 'comment'
        order = 50
        values = [1, "CHARMM force field: No FF information parsed..."]
        return values, format_string, comment, order
    
    def POINTERS(self):
        format_string = '10I8'
        comment = 'comment'
        order = 100
        natom = len(self.topology.atoms)
        ntypes = len(self.topology.atom_types)
        nbonh = len(self.topology.bonds_wH)
        mbona = len(self.topology.bonds_woH)
        ntheth = len(self.topology.angles_wH)
        mtheta = len(self.topology.angles_woH)
        nphih = len(self.topology.dihedrals_wH)
        mphia = len(self.topology.dihedrals_woH)
        nhparm = 0 # unused by amber
        nparm = 0 # 1 if topology file created by addles
        nnb = sum( len(atom.exclusions) if len(atom.exclusions)>0 else 1
                   for atom in self.topology.atoms )
        nres = len( self.topology.residues )
        nbona = mbona 
        ntheta = mtheta
        nphia = mphia
        numbnd = len( self.topology.bond_types )
        numang = len( self.topology.angle_types )
        nptra = len( self.topology.dihedral_types )
        natyp = 0 # unused by amber
        nphb = 0 # for 10-12 H-bond potentials
        ifpert = 0 # unused by current amber
        nbper = 0 # unused by current amber
        ngper = 0 # unused by current amber
        ndper = 0 # unused by current amber
        mbper = 0 # unused by current amber
        mgper = 0 # unused by current amber
        mdper = 0 # unused by current amber
        ifbox = 1 if self.topology.is_periodic else 0
        nmxrs = max( residue.numatoms for residue in self.topology.residues )
        ifcap = 0 
        numextra = 0 # probably zero, check this
        ncopy = 0 # probably zero, check this
    
        values = [
        natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, nhparm,nparm,
        nnb, nres, nbona, ntheta, nphia , numbnd, numang, nptra, natyp, nphb,
        ifpert, nbper, ngper, ndper, mbper, mgper, mdper, ifbox, nmxrs, ifcap,
        numextra , ncopy
        ]
        return values, format_string, comment, order
    
    def ATOM_NAME(self):
        format_string = '20A4'
        comment = 'comment'
        order = 200
        values = [ atom.name for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def CHARGE(self):
        format_string = '3E24.16'
        comment = 'comment'
        order = 300
        k = 18.2223
        values = [ k*atom.charge for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def ATOMIC_NUMBER(self):
        format_string = '10I8'
        comment = NOCOMMENT
        order = 400
        values = [ 1000+i for i,a in enumerate(self.topology.atoms) ]
        return values, format_string, comment, order
    
    def MASS(self):
        format_string = '5E16.8'
        comment = NOCOMMENT
        order = 500
        values = [ atom.mass for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def ATOM_TYPE_INDEX(self):
        format_string = '10I8'
        comment = 'comment'
        order = 600
        values = [ atom.typecode+1 for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def NUMBER_EXCLUDED_ATOMS(self):
        format_string = '10I8'
        comment = NOCOMMENT
        order = 700
        values = [ len(atom.exclusions) if len(atom.exclusions)>0 else 1
                   for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def NONBONDED_PARM_INDEX(self):
        format_string = '10I8'
        comment = 'comment'
        order = 800
        numtypes = len(self.topology.atom_types)
        irange = range(1, numtypes+1)
        values = [ _nb_parm_index(i,j) for i in irange for j in irange ]
        return values, format_string, comment, order
    
    def RESIDUE_LABEL(self):
        format_string = '20A4'
        comment = 'comment'
        order = 900
        values = [ residue.name for residue in self.topology.residues ]
        return values, format_string, comment, order
    
    def RESIDUE_POINTER(self):
        format_string = '10I8'
        comment = 'comment'
        order = 1000
        previous = -1
        values = [ residue.first+1 for residue in self.topology.residues ]
        return values, format_string, comment, order
    
    def BOND_FORCE_CONSTANT(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1100
        values = [ 0.5*bond.k for bond in self.topology.bond_types ]
        return values, format_string, comment, order
    
    def BOND_EQUIL_VALUE(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1200
        values = [ bond.r0 for bond in self.topology.bond_types ]
        return values, format_string, comment, order
    
    def ANGLE_FORCE_CONSTANT(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1300
        values = [ 0.5*angle.k for angle in self.topology.angle_types ]
        return values, format_string, comment, order
    
    def ANGLE_EQUIL_VALUE(self):
        format_string = '3E25.17'
        comment = 'comment'
        order = 1400
        values = [ angle.theta0 for angle in self.topology.angle_types ]
        return values, format_string, comment, order
    
    def DIHEDRAL_FORCE_CONSTANT(self):
        format_string = '5E16.8'
        comment = NOCOMMENT
        order = 1500
        values = [ 0.5*dihedral.k
                    for dihedral in self.topology.dihedral_types ]
        return values, format_string, comment, order
    
    def DIHEDRAL_PERIODICITY(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1600
        values = [ dihedral.n for dihedral in self.topology.dihedral_types ]
        return values, format_string, comment, order
    
    def DIHEDRAL_PHASE(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1700
        values = [ dihedral.phi0
                    for dihedral in self.topology.dihedral_types ]
        return values, format_string, comment, order
    
    def SCEE_SCALE_FACTOR(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1800
        values = [ 1.0 for dihedral in self.topology.dihedral_types ]
        return values, format_string, comment, order
    
    def SCNB_SCALE_FACTOR(self):
        format_string = '5E16.8'
        comment = 'comment'
        order = 1900
        values = [ 1.0 for dihedral in self.topology.dihedral_types ]
        return values, format_string, comment, order
    
    def SOLTY(self):
        format_string = '20A4'
        comment = NOCOMMENT
        order = 2000
        values = []
        return values, format_string, comment, order
    
    def LENNARD_JONES_ACOEF(self):
        format_string = '3E24.16'
        comment = 'comment'
        order = 2100
        values = [ pair.c12 for pair in self.topology.lj_pair_types ]
        return values, format_string, comment, order
    
    def LENNARD_JONES_BCOEF(self):
        format_string = '3E24.16'
        comment = 'comment'
        order = 2200
        values = [ pair.c6 for pair in self.topology.lj_pair_types ]
        return values, format_string, comment, order
    
    def BONDS_INC_HYDROGEN(self):
        format_string = '10I8'
        comment = NOCOMMENT
        order = 2300
        values = _get_amber_indices(self.topology.bonds_wH)
        return values, format_string, comment, order
    
    def BONDS_WITHOUT_HYDROGEN(self):
        format_string = '10I8'
        comment = 'comment'
        order = 2400
        values = _get_amber_indices(self.topology.bonds_woH)
        return values, format_string, comment, order
    
    def ANGLES_INC_HYDROGEN(self):
        format_string = '10I8'
        comment = 'comment'
        order = 2500
        values = _get_amber_indices(self.topology.angles_wH)
        return values, format_string, comment, order
    
    def ANGLES_WITHOUT_HYDROGEN(self):
        format_string = '10I8'
        comment = 'comment'
        order = 2600
        values = _get_amber_indices(self.topology.angles_woH)
        return values, format_string, comment, order
    
    def DIHEDRALS_INC_HYDROGEN(self):
        format_string = '10I8'
        comment = 'comment'
        order = 2700
        values = _get_amber_indices(self.topology.dihedrals_wH)
        return values, format_string, comment, order
    
    def DIHEDRALS_WITHOUT_HYDROGEN(self):
        format_string = '10I8'
        comment = 'comment'
        order = 2800
        values = _get_amber_indices(self.topology.dihedrals_woH)
        return values, format_string, comment, order
    
    def EXCLUDED_ATOMS_LIST(self):
        format_string = '10I8'
        comment = 'comment'
        order = 2900
        values = []
        for atom in self.topology.atoms:
            numexcl = len(atom.exclusions)
            exclusions = [ e+1 for e in atom.exclusions ] if numexcl>0 else [0]
            values.extend(exclusions)
        return values, format_string, comment, order
    
    def HBOND_ACOEF(self):
        format_string = '20A4'
        comment = 'comment'
        order = 3000
        values = []
        return values, format_string, comment, order
    
    def HBOND_BCOEF(self):
        format_string = '20A4'
        comment = 'comment'
        order = 3200
        values = []
        return values, format_string, comment, order
    
    def HBCUT(self):
        format_string = '20A4'
        comment = 'comment'
        order = 3300
        values = []
        return values, format_string, comment, order
    
    def AMBER_ATOM_TYPE(self):
        format_string = '20A4'
        comment = 'comment'
        order = 3400
        atomtypes = self.topology.atom_types
        values = [ atomtypes[atom.typecode] for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def TREE_CHAIN_CLASSIFICATION(self):
        format_string = '20A4'
        comment = 'All items BLA in Chamber topology'
        order = 3500
        values = [ 'BLA' for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def JOIN_ARRAY(self):
        format_string = '10I8'
        comment = 'comment'
        order = 3600
        values = [ 0 for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def IROTAT(self):
        format_string = '10I8'
        comment = 'comment'
        order = 3700
        values = [ 0 for atom in self.topology.atoms ]
        return values, format_string, comment, order
    
    def SOLVENT_POINTERS(self):
        format_string = '3I8'
        comment = 'comment'
        order = 3730
        values = [ self.topology.num_solute_residues,
                   len(self.topology.atoms_per_molecule),
                   self.topology.num_solute_molecules+1
                   ]
        return values, format_string, comment, order

    def ATOMS_PER_MOLECULE(self):
        format_string = '10I8'
        comment = 'comment'
        order = 3730
        values = [ numatoms for numatoms in self.topology.atoms_per_molecule ]
        return values, format_string, comment, order
    
    def CHARMM_UREY_BRADLEY_COUNT(self):
        format_string = '2I8'
        comment = NOCOMMENT
        order = 3800
        values = [0,0]
        return values, format_string, comment, order
    
    def CHARMM_UREY_BRADLEY(self):
        format_string = '10I8'
        comment = NOTERMS
        order = 3900
        values = []
        return values, format_string, comment, order
    
    def CHARMM_UREY_BRADLEY_FORCE_CONSTANT(self):
        format_string = '20A4'
        comment = 'comment'
        order = 4000
        values = []
        return values, format_string, comment, order
    
    def CHARMM_UREY_BRADLEY_EQUIL_VALUE(self):
        format_string = '20A4'
        comment = 'comment'
        order = 4100
        values = []
        return values, format_string, comment, order
    
    def CHARMM_NUM_IMPROPERS(self):
        format_string = 'I8'
        comment = NOCOMMENT
        order = 4200
        values = [ len(self.topology.impropers_woH)
                    + len(self.topology.impropers_wH), ]
        return values, format_string, comment, order
    
    def CHARMM_IMPROPERS(self):
        format_string = "10I8"
        comment = 'comment'
        order = 4300
        impropers = list(self.topology.impropers_wH)
        impropers.extend(self.topology.impropers_woH)
        values = _get_amber_indices(impropers, impropers = True)
        return values, format_string, comment, order
    
    def CHARMM_NUM_IMPR_TYPES(self):
        format_string = 'I8'
        comment = NOCOMMENT
        order = 4400
        values = [ len( self.topology.improper_types ), ]
        return values, format_string, comment, order
    
    def CHARMM_IMPROPER_FORCE_CONSTANT(self):
        format_string = '5E16.8'
        comment = NOCOMMENT
        order = 4500
        values = [ improper.k for improper in self.topology.improper_types ]
        return values, format_string, comment, order
    
    def CHARMM_IMPROPER_PHASE(self):
        format_string = '5E16.8'
        comment = 'In degrees'
        order = 4600
        values = [ improper.xi0 for improper in self.topology.improper_types ]
        return values, format_string, comment, order
    
    def LENNARD_JONES_14_ACOEF(self):
        format_string = '3E24.16'
        comment = 'comment'
        order = 4700
        values = [ pair.c12_14 for pair in self.topology.lj_pair_types ]
        return values, format_string, comment, order
    
    def LENNARD_JONES_14_BCOEF(self):
        format_string = '3E24.16'
        comment = 'comment'
        order = 4800
        values = [ pair.c6_14 for pair in self.topology.lj_pair_types ]
        return values, format_string, comment, order
     
# for bonds, angles, and dihedrals, but not chamber impropers
def _amber_index(index):
    return 3*(index-1) if index>0 else -3*(-index-1)

def _get_amber_indices(interactions, impropers = False):
    values = []
    for interaction in interactions:
        atoms = interaction.atoms
        if interaction.is_excluding_14():
            atoms = list(atoms) # avoid mutating original
            atoms[2] *= -1
        for atom_index in atoms:
            index = atom_index+1
            index = index if impropers else _amber_index(index)
            values.append(index)
        values.append(interaction.typecode+1)
    return values

def _nb_parm_index(i, j):
    #NOTE :
    # while documentation suggests that these indices are up to the user,
    # in fact if using charmm 1-4 parameters this index patterns is
    # hard-coded in ene.F90
    if not type(i): throw(Exception("i must be int"))
    if not type(j): throw(Exception("j must be int"))
    i,j = (i,j) if i<=j else (j,i) # enforce i <= j
    return j*(j-1)/2 + i

def _section(title, comment, format_string, values):
    header = _section_header(title, comment, format_string)
    body = _section_body(values, format_string)
    return header+body

def _section_header(title, comment, format_string):
    flag = "%FLAG {}\n".format(title)
    nocomment = comment == NOCOMMENT
    nocomment = True
    com = "" if nocomment else "%COMMENT {}\n".format(title)
    fmt = "%FORMAT({})\n".format(format_string)
    return flag + com + fmt

def _section_body(values, format_string):
    return FortranWriter(format_string).write(values) + '\n'

