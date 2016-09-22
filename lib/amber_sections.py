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
# amber_sections.py to sections in the topology file, it is important
# that the functions in amber_sections all correspond to sections. Helper
# functions have been put in amber_helpers.py. Unfortunately, we cannot simply
# import all the amber_helpers functions, as inspect will then treat them as
# though they are part of amber_sections.py.

import amber_helpers as h
import output_helpers as oh

# Constants
NOCOMMENT = "" # comment when no comment required
NOTERMS = "This section intentionally left empty."
LINEWIDTH = 80

def CTITLE(data, ifp, united):
    format_string = '20a4'
    comment = NOCOMMENT
    order = 0
    title = data.var["rnme"]
    title = title + " "*(LINEWIDTH - len(title))
    values = [title[i:i+4] for i in range(0, LINEWIDTH, 4)]
    return values, format_string, comment, order

def POINTERS(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 10
    natom = len(h.get_atoms(data, united))
    ntypes = len(set(h.get_atom_types(data, united)))
    nbonh = len(h.filter_hydrogen(h.get_bonds, True, data, united))
    mbona = len(h.filter_hydrogen(h.get_bonds, False, data, united))
    ntheth = len(h.filter_hydrogen(h.get_angles, True, data, united))
    mtheta = len(h.filter_hydrogen(h.get_angles, False, data, united))
    nphih = len(h.filter_hydrogen(h.get_dihedrals, True, data, united))
    mphia = len(h.filter_hydrogen(h.get_dihedrals, False, data, united))
    nhparm = 0 # unused by amber
    nparm = 0 # "used to determine if addles created prmtop" whatever that means
    nnb = sum( len(e) for e in h.get_exclusions(data, united) )
    nres = len( h.get_residue_set(data, united) )
    nbona = mbona 
    ntheta = mtheta
    nphia = mphia
    numbnd = len( h.get_interaction_types(h.get_bonds, data, united) )
    numang = len( h.get_interaction_types(h.get_angles, data, united) )
    nptra = len( h.get_interaction_types(h.get_dihedrals, data, united) )
    natyp = 0 # unused by amber
    nphb = 0 # for 10-12 H-bond potentials
    ifpert = 0 # unused by current amber
    nbper = 0 # unused by current amber
    ngper = 0 # unused by current amber
    ndper = 0 # unused by current amber
    mbper = 0 # unused by current amber
    mgper = 0 # unused by current amber
    mdper = 0 # unused by current amber
    ifbox = 0 # no periodic box by default
    residues = h.get_residues(data, united)
    residue_count = { r : 0 for r in h.get_residue_set(data, united) }
    for r in residues:
        residue_count[r] += 1
    nmxrs = max(residue_count.values())
    ifcap = 0 # probably zero, check this
    numextra = 0 # probably zero, check this
    ncopy = 0 # probably zero, check this

    values = [
    natom, ntypes, nbonh, mbona, ntheth, mtheta, nphih, mphia, nhparm, nparm,
    nnb, nres, nbona, ntheta, nphia , numbnd, numang, nptra, natyp, nphb,
    ifpert, nbper, ngper, ndper, mbper, mgper, mdper, ifbox, nmxrs, ifcap,
    numextra , ncopy
    ]
    return values, format_string, comment, order

def ATOM_NAME(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 20
    values = [atom['symbol'] for atom in h.get_atoms(data, united)]
    return values, format_string, comment, order

def CHARGE(data, ifp, united):
    format_string = '3E24.16'
    comment = 'comment'
    order = 30
    k = 18.2223
    values = [ k*h.get_charge(atom, united) for atom in h.get_atoms(data,united) ]
    return values, format_string, comment, order

def ATOMIC_NUMBER(data, ifp, united):
    format_string = '10i8'
    comment = NOCOMMENT
    order = 40
    values = [ h.get_atomic_number(atom,united)
               for atom in h.get_atoms(data,united) ]
    return values, format_string, comment, order

def MASS(data, ifp, united):
    format_string = '5E16.8'
    comment = NOCOMMENT
    order = 50
    values = [ h.get_mass(atom,united) for atom in h.get_atoms(data,united) ]
    return values, format_string, comment, order

def ATOM_TYPE_INDEX(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 60
    types = h.get_atom_types(data, united)
    types_set = h.get_atom_type_set(data, united)
    values = [ types_set.index(t)+1 for t in types ]
    return values, format_string, comment, order

def NUMBER_EXCLUDED_ATOMS(data, ifp, united):
    format_string = '10i8'
    comment = NOCOMMENT
    order = 70
    values = [ len(e) for e in h.get_exclusions(data,united) ]
    return values, format_string, comment, order

def NONBONDED_PARM_INDEX(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 80
    types_set = h.get_atom_type_set(data, united)
    irange = range(1,len(types_set)+1)
    values = [ h.nb_parm_index(i,j) for i in irange for j in irange ]
    return values, format_string, comment, order

def RESIDUE_LABEL(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 100
    values = h.get_residue_set(data, united)
    return values, format_string, comment, order

def RESIDUE_POINTER(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 120
    residues = h.get_residues(data, united)
    residues_set = h.get_residue_set(data, united)
    values = [ residues.index(r)+1 for r in residues_set ]
    return values, format_string, comment, order

def BOND_FORCE_CONSTANT(data, ifp, united):
    format_string = '5E16.8'
    comment = 'comment'
    order = 130
    types = h.get_interaction_types(h.get_bonds, data, united)
    values = [ types[t]['hfc']*0.5*h.KILOJOULE/h.NANOMETRE**2
               for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def BOND_EQUIL_VALUE(data, ifp, united):
    format_string = '5E16.8'
    comment = 'comment'
    order = 140
    types = h.get_interaction_types(h.get_bonds, data, united)
    values = [ types[t]['value']*h.NANOMETRE for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def ANGLE_FORCE_CONSTANT(data, ifp, united):
    format_string = '5E16.8'
    comment = 'comment'
    order = 150
    types = h.get_interaction_types(h.get_angles, data, united)
    values = [ types[t]['hfc']*0.5*h.KILOJOULE/h.DEGREE**2
               for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def ANGLE_EQUIL_VALUE(data, ifp, united):
    format_string = '3E25.17'
    comment = 'comment'
    order = 160
    types = h.get_interaction_types(h.get_angles, data, united)
    values = [ types[t]['value']*h.DEGREE for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def DIHEDRAL_FORCE_CONSTANT(data, ifp, united):
    format_string = '5E16.8'
    comment = NOCOMMENT
    order = 170
    types = h.get_dihedral_types(data, united)
    values = [ types[t][0]['fc']*0.5*h.KILOJOULE for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def DIHEDRAL_PERIODICITY(data, ifp, united):
    format_string = '5E16.8'
    comment = 'comment'
    order = 180
    types = h.get_dihedral_types(data, united)
    values = [ types[t][0]['mul'] for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def DIHEDRAL_PHASE(data, ifp, united):
    format_string = '5E16.8'
    comment = 'comment'
    order = 190
    types = h.get_dihedral_types(data, united)
    values = [ h.get_phase(types[t][0])*h.DEGREE for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def SCEE_SCALE_FACTOR(data, ifp, united):
    format_string = '5e16.8'
    comment = 'comment'
    order = 200
    dihedral_types = h.get_dihedral_types(data, united)
    values = [ 0.0 if dihedral_types[key][1] else 1.0
               for key in sorted(dihedral_types.keys()) ]
    return values, format_string, comment, order

def SCNB_SCALE_FACTOR(data, ifp, united):
    format_string = '5e16.8'
    comment = 'comment'
    order = 210
    dihedral_types = h.get_dihedral_types(data, united)
    values = [ 0.0 if dihedral_types[key][1] else 1.0
               for key in sorted(dihedral_types.keys()) ]
    return values, format_string, comment, order

def SOLTY(data, ifp, united):
    format_string = '20a4'
    comment = NOCOMMENT
    order = 220
    values = []
    return values, format_string, comment, order

def LENNARD_JONES_ACOEF(data, ifp, united):
    format_string = '3E24.16'
    comment = 'comment'
    order = 230
    values = h.lennard_jones_coef(data, ifp, united, 'A', False)
    return values, format_string, comment, order

def LENNARD_JONES_BCOEF(data, ifp, united):
    format_string = '3E24.16'
    comment = 'comment'
    order = 240
    values = h.lennard_jones_coef(data, ifp, united, 'B', False)
    return values, format_string, comment, order

def BONDS_INC_HYDROGEN(data, ifp, united):
    format_string = '10I8'
    comment = NOCOMMENT
    order = 250
    bonds = h.filter_hydrogen(h.get_bonds, True, data, united)
    values = []
    for (i,bond) in bonds:
        ai,aj = h.get_interaction_atom_indices(bond, data, united)
        values.append(h.amber_index(ai))
        values.append(h.amber_index(aj))
        values.append(i)
    return values, format_string, comment, order

def BONDS_WITHOUT_HYDROGEN(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 260
    bonds = h.filter_hydrogen(h.get_bonds, False, data, united)
    values = []
    for (i,bond) in bonds:
        ai,aj = h.get_interaction_atom_indices(bond, data, united)
        values.append(h.amber_index(ai))
        values.append(h.amber_index(aj))
        values.append(i)
    return values, format_string, comment, order

def ANGLES_INC_HYDROGEN(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 170
    angles = h.filter_hydrogen(h.get_angles, True, data, united)
    values = []
    for (i,angle) in angles:
        ai,aj,ak = h.get_interaction_atom_indices(angle, data, united)
        values.append(h.amber_index(ai))
        values.append(h.amber_index(aj))
        values.append(h.amber_index(ak))
        values.append(i)
    return values, format_string, comment, order

def ANGLES_WITHOUT_HYDROGEN(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 280
    angles = h.filter_hydrogen(h.get_angles, False, data, united)
    values = []
    for (i,angle) in angles:
        ai,aj,ak = h.get_interaction_atom_indices(angle, data, united)
        values.append(h.amber_index(ai))
        values.append(h.amber_index(aj))
        values.append(h.amber_index(ak))
        values.append(i)
    return values, format_string, comment, order

def DIHEDRALS_INC_HYDROGEN(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 290
    dihedrals = h.filter_hydrogen(h.get_dihedrals, True, data, united)
    values = []
    for (i,dihedral) in dihedrals:
        ai,aj,ak, al = h.get_interaction_atom_indices(dihedral, data, united)
        values.append(h.amber_index(ai))
        values.append(h.amber_index(aj))
        values.append(h.amber_index(ak))
        values.append(h.amber_index(al))
        values.append(i)
    return values, format_string, comment, order

def DIHEDRALS_WITHOUT_HYDROGEN(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 300
    dihedrals = h.filter_hydrogen(h.get_dihedrals, False, data, united)
    values = []
    for (i,dihedral) in dihedrals:
        ai,aj,ak, al = h.get_interaction_atom_indices(dihedral, data, united)
        values.append(h.amber_index(ai))
        values.append(h.amber_index(aj))
        values.append(h.amber_index(ak))
        values.append(h.amber_index(al))
        values.append(i)
    return values, format_string, comment, order

def EXCLUDED_ATOMS_LIST(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 305
    values = []
    [ values.extend(e) for e in h.get_exclusions(data, united) ]
    return values, format_string, comment, order

def HBOND_ACOEF(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 310
    values = []
    return values, format_string, comment, order

def HBOND_BCOEF(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 320
    values = []
    return values, format_string, comment, order

def HBCUT(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 330
    values = []
    return values, format_string, comment, order

def AMBER_ATOM_TYPE(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 340
    values = h.get_atom_type_set(data, united)
    return values, format_string, comment, order

def TREE_CHAIN_CLASSIFICATION(data, ifp, united):
    format_string = '20a4'
    comment = 'All items BLA in Chamber topology'
    order = 350
    values = [ 'BLA' for atom in h.get_atoms(data, united) ]
    return values, format_string, comment, order

def JOIN_ARRAY(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 351
    values = [0]*len(h.get_atoms(data, united))
    return values, format_string, comment, order

def IROTAT(data, ifp, united):
    format_string = '10I8'
    comment = 'comment'
    order = 352
    values = [0]*len(h.get_atoms(data, united))
    return values, format_string, comment, order

def FORCE_FIELD_TYPE(data, ifp, united):
    format_string ='i2,a78'
    comment = 'comment'
    order = 15 #420
    values = [1, "CHARMM force field: No FF information parsed..."]
    return values, format_string, comment, order

def CHARMM_UREY_BRADLEY_COUNT(data, ifp, united):
    format_string = '2i8'
    comment = NOCOMMENT
    order = 430
    values = [0,0]
    return values, format_string, comment, order

def CHARMM_UREY_BRADLEY(data, ifp, united):
    format_string = '10I8'
    comment = NOTERMS
    order = 440
    values = []
    return values, format_string, comment, order

def CHARMM_UREY_BRADLEY_FORCE_CONSTANT(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 450
    values = []
    return values, format_string, comment, order

def CHARMM_UREY_BRADLEY_EQUIL_VALUE(data, ifp, united):
    format_string = '20a4'
    comment = 'comment'
    order = 460
    values = []
    return values, format_string, comment, order

def CHARMM_NUM_IMPROPERS(data, ifp, united):
    format_string = 'i8'
    comment = NOCOMMENT
    order = 470
    values = [ len(h.get_impropers(data, united)), ]
    return values, format_string, comment, order

def CHARMM_IMPROPERS(data, ifp, united):
    format_string = "10i8"
    comment = 'comment'
    order = 480
    impropers = h.get_impropers(data, united)
    types =sorted(h.get_improper_types(data,united).keys())
    type_index = { t:i+1 for i,t in enumerate(types) }
    values = []
    for (i,improper) in enumerate(h.get_impropers(data, united)):
        ai,aj,ak,al = h.get_interaction_atom_indices(improper, data, united)
        # NOTE DIFFERENCE HERE COMPARED TO DIHEDRALS!
        # here we use the normal indices, not that calculated by amber_index(). 
        # It seems the CHAMBER devs like to be different.
        values.append(ai)
        values.append(aj)
        values.append(ak)
        values.append(al)
        values.append(type_index[improper['code']])
    return values, format_string, comment, order

def CHARMM_NUM_IMPR_TYPES(data, ifp, united):
    format_string = 'i8'
    comment = NOCOMMENT
    order = 490
    values = [ len( h.get_improper_types(data, united) ), ]
    return values, format_string, comment, order

def CHARMM_IMPROPER_FORCE_CONSTANT(data, ifp, united):
    format_string = '5E16.8'
    comment = NOCOMMENT
    order = 500
    types = h.get_improper_types(data, united)
    values = [ types[t]['fc']*0.5*h.KILOJOULE/h.DEGREE**2
               for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def CHARMM_IMPROPER_PHASE(data, ifp, united):
    format_string = '5E16.8'
    comment = 'In degrees'
    order = 510
    types = h.get_improper_types(data, united)
    print(types)
    values = [ types[t]['value']
               for t in sorted(types.keys()) ]
    return values, format_string, comment, order

def LENNARD_JONES_14_ACOEF(data, ifp, united):
    format_string = '3E24.16'
    comment = 'comment'
    order = 513
    one4 = True
    values = h.lennard_jones_coef(data, ifp, united, 'A', one4)
    return values, format_string, comment, order

def LENNARD_JONES_14_BCOEF(data, ifp, united):
    format_string = '3E24.16'
    comment = 'comment'
    order = 518
    one4 = True
    values = h.lennard_jones_coef(data, ifp, united, 'B', one4)
    return values, format_string, comment, order
 
