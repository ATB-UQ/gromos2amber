from .GromosTopologyParser import GromosTopologyParser
from math import sqrt

KILOJOULE = 1.0/4.184 # kCal
NANOMETRE = 10.0 # angstroms
DEGREE = 3.141592653589793/180.0 #radians

class Topology:

    # _wH := with hydrogen
    # _woH := without hydrogen

    def __init__(self, io):
        gromos = GromosTopologyParser(io)

        self.title = gromos.TITLE()

        self.atoms, self.residues  = _read_atoms_and_residues(gromos)
        self.atom_types = _read_atom_types(gromos)

        self.lj_pair_types = _read_lj_pair_types(gromos)

        self.bond_types = _read_bond_types(gromos)
        self.angle_types = _read_angle_types(gromos)
        self.dihedral_types = _read_dihedral_types(gromos)
        self.improper_types = _read_improper_types(gromos)

        self.dihedral_types.append(DihedralType(0.0,0.0,1.0)) #dummy for 1-4

        ri = _read_bonded_interaction
        self.bonds_woH     = ri(gromos.BOND(H=False))
        self.angles_woH    = ri(gromos.BONDANGLE(H=False))
        self.dihedrals_woH = ri(gromos.DIHEDRAL(H=False))
        self.impropers_woH = ri(gromos.IMPDIHEDRAL(H=False))

        self.bonds_wH     = ri(gromos.BOND(H=True))
        self.angles_wH    = ri(gromos.BONDANGLE(H=True))
        self.dihedrals_wH = ri(gromos.DIHEDRAL(H=True))
        self.impropers_wH = ri(gromos.IMPDIHEDRAL(H=True))
        # Marks dihedrals for which 1-4 interactions must be excluded
        all_dihedrals = self.dihedrals_wH + self.dihedrals_woH
        _fix_14_exclusions(self.atoms, all_dihedrals)
        # Add dummy dihedrals to force 1-4 interactions where required
        self.dihedrals_woH.extend(_extra_dihedrals(
            self.atoms,
            self.dihedrals_wH,
            self.dihedrals_woH,
            len(self.dihedral_types)-1),
        )
        self.is_periodic = True
        one_on_4_pi_eps0 = gromos.PHYSICALCONSTANTS()[0] *KILOJOULE*NANOMETRE
        self.charge_prefactor = sqrt(one_on_4_pi_eps0) #gromos

        self.atoms_per_molecule = _read_atoms_per_solute_molecule(gromos)
        self.num_solute_molecules = len(self.atoms_per_molecule)

        self.solvent_atoms = _read_solvent(gromos)
        bondinfo = _read_solvent_bonds(gromos, len(self.bond_types))
        self.solvent_bonds, self.solvent_bond_types = bondinfo

    def add_solvent(self, num_solvent_molecules, residue_name):

        atoms_per_solvent = len(self.solvent_atoms)
        num_solute_atoms = len(self.atoms)

        self.atoms_per_molecule.extend([
            len(self.solvent_atoms)
            for i in range(num_solvent_molecules)
        ])

        self.num_solute_residues = len(self.residues)

        self.residues.extend([
            Residue(residue_name,
                atoms_per_solvent*i + num_solute_atoms,
                atoms_per_solvent)
            for i in range(num_solvent_molecules)
        ])

        self.bond_types.extend(self.solvent_bond_types)

        all_solvent_atoms, all_solvent_bonds = _create_solvent(
            self.solvent_atoms,
            self.solvent_bonds,
            num_solvent_molecules,
            num_solute_atoms
        )

        self.atoms.extend(all_solvent_atoms)
        self.bonds_wH.extend(all_solvent_bonds)


    def get_title(self): return self.title.replace('\n','_')

##### End of Topology class #####

def _read_solvent_bonds(gromos, num_bond_types):
    ii,jj,lengths = gromos.SOLVENTCONSTR()
    lengths = [ r0*NANOMETRE for r0 in lengths ]
    solvent_bond_types = [ BondType(0.0,r0) for r0 in set(lengths) ]
    typecodes = {
        bondtype.r0 : index+num_bond_types
        for index,bondtype in enumerate(solvent_bond_types)
    }
    solvent_bonds = [
        Interaction([i-1,j-1], typecodes[r0])
        for i,j,r0 in zip(ii,jj,lengths)
    ]
    
    return solvent_bonds, solvent_bond_types



def _read_atoms_and_residues(gromos):
    data = gromos.SOLUTEATOM()
    atomindex, residue_number, name, typecode = data[0:4]
    mass, charge, _, exclusions, neigh14 = data[4: ]
    numatoms = len(atomindex)
    atoms = [
        Atom(
            name[i],
            typecode[i]-1, 
            mass[i],
            charge[i],
            [ e-1 for e in exclusions[i] ],
            [ n14-1 for n14 in neigh14[i] ],
        )
        for i in range(numatoms)
    ]
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
    return [
        name[:4] if len(name.strip())>4 else name.strip()
        for i,name in enumerate(names)
    ]

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
    return [
        DihedralType(k*unitk,phi0*unitphi0,n)
        for k,phi0,n in zip(spring, phase, multiplicity)
    ]

def _read_improper_types(gromos):
    spring,angle = gromos.IMPDIHEDRALTYPE()
    unitk = KILOJOULE/DEGREE**2
    unitxi0 = DEGREE
    return [ ImproperType(k*unitk,xi0*unitxi0)
                for k,xi0 in zip(spring, angle) ]

# processes bond, angle, dihedral, improper columns from gromos parser
def _read_bonded_interaction(gromos_columns):
    cols = gromos_columns
    numcols = len(gromos_columns) - 1
    numrows = len(gromos_columns[0])
    return [
        Interaction([ cols[c][r]-1 for c in range(numcols) ], cols[-1][r]-1 )
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
    return [
        LJPairType(
            typei[i]-1,
            typej[i]-1,
            c12[i]*unit12,
            c6[i]*unit6,
            c12_14[i]*unit12,
            c6_14[i]*unit6,
        )
        for i in range(numpairs)
    ]

def _read_solvent(gromos):
    _,name,typecode,mass,charge = gromos.SOLVENTATOM()
    n = len(name)
    atoms = [ Atom(name[i].strip(),
                   typecode[i]-1, 
                   mass[i],
                   charge[i],
                   [], [])
              for i in range(n) ]
    return atoms

def _create_solvent(solvent_atoms, solvent_bonds,
                    num_molecules, first_solvent_index):
    name = [ atom.name for atom in solvent_atoms ]
    m = [ atom.mass for atom in solvent_atoms ]
    q = [ atom.charge for atom in solvent_atoms ]
    tc = [ atom.typecode for atom in solvent_atoms ]
    n = len(solvent_atoms)
    i0 = first_solvent_index
    atoms = [
        Atom(name[i%n],
            tc[i%n], 
            m[i%n],
            q[i%n],
            [ j+i0 for j in range(int(i/n)*n,int(i/n+1)*n) if not j==i ],
            [],
        )
        for i in range(num_molecules * n)
    ]

    bonds = []
    for mol in range(num_molecules):
        for bond in solvent_bonds:
            bi = i0+mol*n+bond.atoms[0]
            bj = i0+mol*n+bond.atoms[1]
            bonds.append(Interaction([bi,bj],bond.typecode))

    return atoms, bonds

def _read_atoms_per_solute_molecule(gromos):
    mol_last_index = gromos.SOLUTEMOLECULES()
    nummol = len(mol_last_index)
    return [
        mol_last_index[i] - (0 if i==0 else mol_last_index[i-1])
        for i in range(nummol)
    ]

def _extra_dihedrals(atoms, dihedrals_wH, dihedrals_woH, dummy_typecode):
    extra = []
    all_dihedrals = list(dihedrals_wH)
    all_dihedrals.extend(dihedrals_woH)
    found = [ { n14:False for n14 in atom.neigh14 } for atom in atoms ]
    for dihedral in all_dihedrals:
        i, l = dihedral.atoms[0], dihedral.atoms[3]
        i, l = (i,l) if i<l else (l,i)
        found[i][l] = True
    for i,atom in enumerate(atoms):
        for l in atom.neigh14:
           if not found[i][l]:
               extra.append(Interaction([i,i,l,l],dummy_typecode))
    return extra

def _fix_14_exclusions(atoms, dihedrals):
    prev_l = [ [] for atom in atoms ]
    for d,dihedral in enumerate(dihedrals):
        i, l = dihedral.atoms[0], dihedral.atoms[3]
        i, l = (i,l) if i<l else (l,i)
        if l in atoms[i].exclusions_wo14 or l in prev_l[i]:
            dihedral.exclude_14()
        prev_l[i].append(l)

class Atom:
    def __init__(self, name, typecode, mass, charge, exclusions, neigh14):
        self.name = name
        self.typecode = typecode
        self.mass = mass
        self.charge = charge
        self.neigh14 = neigh14
        self.exclusions_wo14 = exclusions
        all_exclusions = list(exclusions)
        all_exclusions.extend(neigh14)
        self.exclusions = sorted(set(all_exclusions))

class Interaction:
    """ Bond, angle, proper dihedral, or improper dihedral """
    def __init__(self, atoms, typecode):
        self.atoms = atoms
        self.typecode = typecode
        self._exclude14 = False
        if len(atoms)==4 and (atoms[2] == 0 or atoms[3] == 0):
            self.atoms.reverse()

    def exclude_14(self, exclude = True):
        if len(self.atoms) != 4: raise Exception("not a dihedral")
        self._exclude14 = exclude

    def is_excluding_14(self): return self._exclude14

class BondType:
    def __init__(self, k, r0):
        self.k, self.r0 = k, r0

class LJPairType:
    def __init__(self, itype, jtype, c12, c6, c12_14, c6_14):
        if c12_14 != 0 and c12 == 0:
            c12 = 1.0e-99
        self.itype, self.jtype = itype, jtype
        self.c12, self.c6, self.c12_14, self.c6_14 = c12, c6, c12_14, c6_14

class AngleType:
    def __init__(self, k, theta0):
        self.k, self.theta0 = k, theta0

class DihedralType:
    def __init__(self, k, phi0, n):
        self.k, self.phi0, self.n = k, phi0, n

class ImproperType:
    def __init__(self, k, xi0):
        self.k, self.xi0 = k, xi0

class Residue:
    def __init__(self, name, first, numatoms):
        self.name, self.first, self.numatoms = name, first, numatoms
