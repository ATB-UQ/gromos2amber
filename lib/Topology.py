from GromosTopologyParser import GromosTopologyParser

KILOJOULE = 1.0/4.184 # kCal
NANOMETRE = 10.0 # angstroms
DEGREE = 3.141592653589793/180.0 #radians

class Topology:

    # _wH := with hydrogen
    # _woH := without hydrogen

    def __init__(self, io, configuration):
        gromos = GromosTopologyParser(io)
        numatoms = len(configuration.positions)

        self.title = gromos.TITLE()

        self.atoms, self.residues  = _read_atoms_and_residues(gromos)
        self.atom_types = _read_atom_types(gromos)

        self.lj_pair_types = _read_lj_pair_types(gromos)

        self.bond_types = _read_bond_types(gromos)
        self.angle_types = _read_angle_types(gromos)
        self.dihedral_types = _read_dihedral_types(gromos)
        self.improper_types = _read_improper_types(gromos)

        ri = _read_bonded_interaction
        self.bonds_woH     = ri(gromos.BOND(H=False))
        self.angles_woH    = ri(gromos.BONDANGLE(H=False))
        self.dihedrals_woH = ri(gromos.DIHEDRAL(H=False))
        self.impropers_woH = ri(gromos.IMPDIHEDRAL(H=False))

        self.bonds_wH     = ri(gromos.BOND(H=True))
        self.angles_wH    = ri(gromos.BONDANGLE(H=True))
        self.dihedrals_wH = ri(gromos.DIHEDRAL(H=True))
        self.impropers_wH = ri(gromos.IMPDIHEDRAL(H=True))

        self.is_periodic = True

        numsolvent = numatoms - len(self.atoms)
        solvent_atoms, atomspersolvent = _read_solvent(gromos, numsolvent)

    def get_title(self): return self.title.replace('\n','_')


##### End of Topology class #####


def _read_atoms_and_residues(gromos):
    data = gromos.SOLUTEATOM()
    atomindex, residue_number, name, typecode = data[0:4]
    mass, charge, _, exclusions, neigh14 = data[4: ]
    numatoms = len(atomindex)
    atoms = [ Atom(name[i],
                        typecode[i]-1, 
                        mass[i],
                        charge[i],
                        exclusions[i],
                        neigh14[i])
                    for i in range(numatoms) ]
    residue_names = gromos.RESNAMES()
    residues = []
    previous = -1
    for i,res in enumerate(residue_number):
        if res != previous:
            residues.append( Residue(residue_names[res-1], i, 1) )
            previous = res
        else:
            residues[-1].numatoms += 1
    return atoms, residues


def _read_atom_types(gromos):
    names = gromos.ATOMTYPENAME()
    # Amber limits type names to 2 characters.
    # Use typecode instead when this is violated
    return [ str(i+1) if len(name.strip())>2 else name.strip()
             for i,name in enumerate(names) ]

def _read_bond_types(gromos):
    _,spring,length = gromos.BONDSTRETCHTYPE()
    unitk = KILOJOULE/NANOMETRE**2
    unitr0 = NANOMETRE
    return [ BondType(k*unitk,r0*unitr0) for k,r0 in zip(spring, length) ]

def _read_angle_types(gromos):
    _,spring,angle = gromos.BONDANGLEBENDTYPE()
    unitk = KILOJOULE/DEGREE**2
    unitt0 = DEGREE
    return [ AngleType(k*unitk,t0*unitt0) for k,t0 in zip(spring, angle) ]

def _read_dihedral_types(gromos):
    spring, phase, multiplicity  = gromos.TORSDIHEDRALTYPE()
    unitk = KILOJOULE
    unitphi0 = DEGREE
    return [ DihedralType(k*unitk,phi0*unitphi0,n)
                for k,phi0,n in zip(spring, phase, multiplicity) ]

def _read_improper_types(gromos):
    spring,angle = gromos.IMPDIHEDRALTYPE()
    unitk = KILOJOULE/DEGREE**2
    unitxi0 = 1.0
    return [ ImproperType(k*unitk,xi0*unitxi0)
                for k,xi0 in zip(spring, angle) ]

# processes bond, angle, dihedral, improper columns from gromos parser
def _read_bonded_interaction(gromos_columns):
    cols = gromos_columns
    numcols = len(gromos_columns) - 1
    numrows = len(gromos_columns[0])
    return [ Interaction([ cols[c][r]-1 for c in range(numcols) ], cols[-1][r])
            for r in range(numrows)
            ]

def _read_lj_pair_types(gromos):
    typei, typej, c12, c6, c12_14, c6_14 = gromos.LJPARAMETERS()
    numpairs = len(typei)
    for p in range(numpairs):
        if p+1 != int(typej[p]*(typej[p]-1)/2+typei[p]):
            raise Exception("bad pair ordering")
    unit6 = KILOJOULE*NANOMETRE**6
    unit12 = KILOJOULE*NANOMETRE**12
    return [ LJPairType(typei[i]-1, typej[i]-1, c12[i]*unit12, c6[i]*unit6,
                        c12_14[i]*unit12, c6_14[i]*unit6)
                for i in range(numpairs)
                ]

def _read_solvent(gromos, numsolvent):
    _,name,typecode,mass,charge = gromos.SOLVENTATOM()
    numatoms = len(name)
    atoms = [ Atom(name[i%numatoms],
                        typecode[i%numatoms]-1, 
                        mass[i%numatoms],
                        charge[i%numatoms],
                        [],
                        [])
                    for i in range(numsolvent) ]
    return atoms, numatoms



class Atom:
    def __init__(self, name, typecode, mass, charge, exclusions, neigh14):
        self.name = name
        self.typecode = typecode
        self.mass = mass
        self.charge = charge
        self.exclusions = exclusions
        self.neigh14 = neigh14

    def add_exclusion(self, exclusion):
        self.exclusions.append(exclusion)

class Interaction:
    """ Bond, angle, proper dihedral, or improper dihedral """
    def __init__(self, atoms, typ):
        self.atoms = atoms
        self.typ = typ

class BondType:
    def __init__(self, k, r0):
        self.k, self.r0 = k, r0

class LJPairType:
    def __init__(self, itype, jtype, c12, c6, c12_14, c6_14):
        self.itype, self.jtype = itype, jtype
        self.c12, self.c6, self.c12_14, self.c6_14 = c12, c6, c12_14, c6_14

class AngleType:
    def __init__(self, k, theta0):
        self.k, self.theta0 = k, theta0

class DihedralType:
    def __init__(self, k, phi0, n):
        self.k, self.phi0, self.n = k, phi0, n
        self._exclude14 = False

    def exclude_14(self, exclude = True):
        self._exclude14 = exclude

    def is_excluding_14(self): return self._exclude14

class ImproperType:
    def __init__(self, k, xi0):
        self.k, self.xi0 = k, xi0

class Residue:
    def __init__(self, name, first, numatoms):
        self.name, self.first, self.numatoms = name, first, numatoms
