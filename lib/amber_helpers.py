""" Helper functions for amber_sections.py """
from algorithm.atb.helpers.parse_gromos import gen_valid_lines
from math import ceil, sqrt
import algorithm.atb._outputs.output_helpers as helpers

KILOJOULE = 1.0/4.184 # kCal
NANOMETRE = 10.0 # angstroms
DEGREE = 3.141592653589793/180.0 #radians

def get_united_atom_prefix(united): return 'u' if united else ''

def get_property(atom, key, united, fallback = False,
                default = False, default_value = None):
    newkey = get_united_atom_prefix(united) + key
    newkey = key if fallback and not newkey in atom else newkey
    return default_value if default and not newkey in atom else atom[newkey]

def get_atoms(data, united):
    atoms = [ atom for atom in data.atoms.values()
                if (not united or 'uindex' in atom) ] 
    sortkey = 'uindex' if united else 'index'
    atoms.sort(key=lambda x: x[sortkey])
    return atoms

def get_charge(atom, united):
    return get_property(atom, 'charge', united, True, True, 0.0)

def get_mass(atom, united): return get_property(atom, 'mass', united, True)

def get_index(atom, united): return get_property(atom, 'index', united, True)

def get_coord(atom):
    return atom[ 'ocoord' if 'ocoord' in atom.keys() else 'coord' ]

def get_atomic_number(atom, united): return 1

def get_ljsym(atom, united): return get_property(atom, 'ljsym', united, True)

def get_residue(atom): return atom['group']

def get_residues(data, united):
    return [ get_residue(atom) for atom in get_atoms(data, united) ]

def get_residue_set(data, united):
    residues = get_residues(data, united)
    return sorted(set(residues), key=residues.index)  

def get_atom_types(data, united):
    return [ helpers.a_ljsym(atom, get_united_atom_prefix(united))
                for atom in get_atoms(data, united) ]

def get_atom_type_set(data, united):
    types = get_atom_types(data, united)
    return sorted(set(types), key=types.index)  

def get_exclusions(data, united):
    atoms = get_atoms(data, united)
    excl= [ [ j for j in sorted(get_property(atom, 'excl', united))
              if j>i+1 ]
            for i,atom in enumerate(atoms) ]

    # exclude 1-4 proper dihedral interactions
    index1 = [   get_index(data.atoms[ dihedral['atoms'][0] ], united)
                 for dihedral in get_dihedrals(data, united)   ]
    index4 = [   get_index(data.atoms[ dihedral['atoms'][3] ], united)
                 for dihedral in get_dihedrals(data, united)   ]
    for i,j in zip(index1, index4):
        ii,jj = (i,j) if i<j else (j,i)
        if not jj in excl[ii-1]:
            excl[ii-1].append(jj)
            excl[ii-1].sort()

    return [ e if len(e)>0 else [0] for e in excl ]

def nb_parm_index(i, j):
    #NOTE :
    # while documentation suggests that these indices are up to the user,
    # in fact if using charmm 1-4 parameters this index patterns is
    # hard-coded in ene.F90
    if not type(i): throw(Exception("i must be int"))
    if not type(j): throw(Exception("j must be int"))
    i,j = (i,j) if i<=j else (j,i) # enforce i <= j
    return j*(j-1)/2 + i

def _get_ifp_lj_types(ifp):
    lines = [ l for l in gen_valid_lines(ifp['SINGLEATOMLJPAIR']) ]
    numtypes = int(lines[0])
    params = []
    num_matrix_lines = int(ceil(numtypes*1.0/20))
    lines_per_type = num_matrix_lines + 2
    for i in range(numtypes):
        words = lines[1+i*lines_per_type].split()
        words2 = lines[2+i*lines_per_type].split()
        matrix = []
        for j in range(3,3+num_matrix_lines):
            matrix.extend([int(word)
                           for word in lines[j+i*lines_per_type].split() ] )

        p = { "IAC"      :int(words[0]),
              "TYPE"     :words[1],
              "SQRTc6"   :float(words[2]),
              "SQRTc12"  :[ float(words[3]), float(words[4]), float(words[5]) ],
              "SQRT14c6" :float(words2[0]),
              "SQRT14c12":float(words2[1]),
              "MATRIX"   :matrix  }
        params.append(p)

    return params

def get_lj_params(ifp):
    lj = _get_ifp_lj_types(ifp)
    params = {}
    types = [ p["TYPE"] for p in lj ]
    params = { ti:{ tj:{} for tj in types } for ti in types }
    for i,ti in enumerate(types):
        for j,tj in enumerate(types):
            C6ij_14  = lj[i]["SQRT14c6"]  * lj[j]["SQRT14c6"]
            C12ij_14 = lj[i]["SQRT14c12"] * lj[j]["SQRT14c12"]
            C6ij     = lj[i]["SQRTc6"]    * lj[j]["SQRTc6"]
            ii = lj[i]["MATRIX"][j]-1
            jj = lj[j]["MATRIX"][i]-1 
            sqrtC12i = lj[i]["SQRTc12"][  lj[i]["MATRIX"][j]-1  ]
            sqrtC12j = lj[j]["SQRTc12"][  lj[j]["MATRIX"][i]-1  ]
            C12ij = sqrtC12i*sqrtC12j
            params[ti][tj]["C6"]     = C6ij
            params[ti][tj]["C12"]    = C12ij
            params[ti][tj]["C6_14"]  = C6ij_14
            params[ti][tj]["C12_14"] = C12ij_14
    return params 

def lennard_jones_coef(data, ifp, united, coef, onefour):
    if not coef in ('A', 'B'):
        throw(Exception("Bad coef argument: "+coef+". Must be A or B"))
    types = get_atom_type_set(data, united)
    numtypes = len(types)
    nb_params = get_lj_params(ifp)
    key, power = ('C6',6) if coef == 'B' else ('C12',12)
    key += "_14" if onefour else ""
    coefs = [ 0.0 for xx in range(numtypes*(numtypes+1)/2) ]
    for i in range(numtypes):
        for j in range(i,numtypes):
            index = nb_parm_index(i+1,j+1)-1
            coefs[index] = nb_params[types[i]][types[j]][key]
            coefs[index] *= KILOJOULE*NANOMETRE**power
    return coefs

def get_bonds(data, united):
    atom_ids = [ atom['id'] for atom in get_atoms(data, united) ]
    return [ bond for bond in data.bonds
             if set(bond['atoms']).issubset(atom_ids) ]

def get_interaction_types(get_interactions, data, united):
    return { interaction['code'][0]['code'] : interaction['code'][0]
             for interaction in get_interactions(data, united) }

def get_angles(data, united):
    atom_ids = [ atom['id'] for atom in get_atoms(data, united) ]
    return [ angle for angle in data.angles
             if set(angle['atoms']).issubset(atom_ids) ]

def isessential(dihedral):
    return not ('essential' in dihedral and not dihedral['essential'])

def get_dihedrals(data, united):
    atom_ids = [ atom['id'] for atom in get_atoms(data, united) ]
    return [ dih for dih in data.dihedrals
             if isessential(dih) and set(dih['atoms']).issubset(atom_ids) ]

def get_dihedral_types(data, united):
    dihedrals = get_dihedrals(data, united)
    dihedral_types = {}
    for dihedral in dihedrals:
        code = dihedral['code'][0]['code']
        atom1 = dihedral['atoms'][0]
        atom4 = dihedral['atoms'][3]
        atom1, atom4 = (atom1, atom4) if atom1<atom4 else (atom4, atom1)
        excl1 = get_property(data.atoms[atom1], 'excl', united)
        exclude14 = atom4 in excl1
        key = (code, exclude14)
        dihedral_types[key] = (dihedral['code'][0], exclude14)
    return dihedral_types

def get_impropers(data, united):
    atom_ids = [ atom['id'] for atom in get_atoms(data, united) ]
    return [ improper for improper in data.impropers
             if set(improper['atoms']).issubset(atom_ids) ]

def get_improper_types(data, united):
    impropers = get_impropers(data, united)
    codes = [ improper['code'] for improper in impropers ]
    return { code : helpers.GROMOS_IMPROPER_DIHEDRALS[code] for code in codes }

def get_phase(dihedral_type):
    #return dihedral_type['phase'] if 'phase' in dihedral_type else 0.0
    return dihedral_type['value']

def filter_hydrogen_test(get_interactions, withH, data, united):
    """ Interactions could be bonds, angles, etc"""
    interactions = get_interactions(data, united)
    types = get_interaction_types(get_interactions, data, united)
    type_index = { t:i+1 for i,t in enumerate(sorted(types.keys())) }
    
    if withH: return []

    return [ (type_index[interact['code'][0]['code']] ,interact)
             for (i,interact) in enumerate(interactions) ]

def filter_hydrogen(get_interactions, withH, data, united):
    """ Interactions could be bonds, angles, etc"""
    interactions = get_interactions(data, united)
    types = get_interaction_types(get_interactions, data, united)
    type_index = { t:i+1 for i,t in enumerate(sorted(types.keys())) }
    return [ (type_index[interact['code'][0]['code']] ,interact)
             for (i,interact) in enumerate(interactions)
             if withH ^ (not 'H' in [ data.atoms[aid]['type']
                                      for aid in interact['atoms'] ]) ]

def get_interaction_atom_indices(interaction, data, united):
        return [ get_property(data.atoms[atomid], 'index', united)
                 for atomid in interaction['atoms'] ]
             
# for bonds, angles, dihedrals, impropers
def amber_index(index): return 3*(index-1)

