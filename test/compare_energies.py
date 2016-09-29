#!/usr/bin/env python

# Thomas Lee August 2016
# contact@tomlee.com.au

# Prints a table comparing the energy outputs files from a gromacs and
# amber MD run

# Amber energies are read from an mdinfo file and gromacs energies from
# the output of gmxdump eg `gmxdump -e ener.edr > ener.edr.dump`

import sys
if len(sys.argv) != 3:
    print("usage: {} mdinfo ener.edr.dump\n".format(sys.argv[0]))
    exit(1)

# mdinfo file from sander
amber_energy_file = sys.argv[1]
# output of `gmxdump -e ener.edr`
gromos_energy_file = sys.argv[2]

# get amber energy info from file

with open(amber_energy_file) as f:
    input = f.read()
terms = input.replace("1-4 ","1-4_").replace('**************','-1e99').split()

# corresponding energy labels in gromacs and amber files
amber_labels = {
        "Etot"   : "Total" , "EPtot"  : "Potential" ,
        "BOND"   : "Bond"      , "ANGLE"  : "Angle" ,
        "DIHED"  : "Dihedral"  , "IMP"    : "Improper" ,
        "1-4_NB" : "LJ-14"        , "1-4_EEL": "Coulomb-14" ,
        "VDWAALS": "LJ (SR)"      , "EELEC"  : "Coulomb (SR)"
            }
gromos_labels = {
        "Total"   : "Total" , "Potential"  : "Potential" ,
        "Bonds"   : "Bond"      , "Angles"  : "Angle" ,
        "Dihedral"  : "Dihedral"  , "Improper"    : "Improper" ,
        "Vdw-14" : "LJ-14"        , "El-14": "Coulomb-14" ,
        "Vdw": "LJ (SR)"      , "El (RF)"  : "Coulomb (SR)"
            }
amber = {label : -9999 for label in amber_labels.values()}
gromos = {label : -9999 for label in gromos_labels.values()}

calorie = 4.184 # J
for i in range(0, len(terms), 3): # read amber energy data
    if terms[i] in amber_labels.keys():
        label = amber_labels[terms[i]]
        value = float(terms[i+2])
        amber[label] = value*calorie

#amber["Coulomb (SR)"] = amber["Coulomb (SR)"] - amber["Coulomb-14"]
#amber["LJ (SR)"] = amber["LJ (SR)"] - amber["LJ-14"]

num_energy_lines = 38
with open(gromos_energy_file) as f: # read gromacs energy data
    gromoslines = f.readlines()
for l,line in enumerate(gromoslines):
    if line == "ENERGIES\n":
        energylines = gromoslines[l+1:l+num_energy_lines]
        break
for line in energylines:
    key, valuestr = line.split(':')
    key = key[2:].strip()
    value = float(valuestr)
    if key in gromos_labels:
        gromos[gromos_labels[key]] = float(value)


fmt = "{: >15s} {: > 14.5e} {: > 14.5e} {: > 10.5e} {: > 10.6f}"
lfmt = "{:15s} {:14s} {:14s} {:5s}"
print( lfmt.format(" quantity", " gromos", " amber", " (amber-gromos)", " )") )
for label in gromos_labels.values():
    if gromos[label] + amber[label] == 0 :
        if gromos[label] == amber[label]:
            dif = 0.
        else:
            dif = 1.
    else:
        dif = (amber[label]-gromos[label])/gromos[label] if gromos[label] != 0 else 999
    x = amber[label]-gromos[label]

    print(fmt.format(label, gromos[label], amber[label], x, dif))

#print("debug")
#print((-amber["LJ-14"])**(-1.0/6))
