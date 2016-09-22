
import sys, os
from inspect import getmembers, isfunction
import amber_sections
import amber_helpers as helper
from fortranformat import FortranRecordWriter as FortranWriter
from math import sqrt

DEBUG_VS_PDB = False

def parm7(data, ifp, united, io):
    section_functions = getmembers(amber_sections, predicate=isfunction)
    sections = []
    ordering = []

    for title,func in section_functions:
        print(title)
        values, format_string, comment, order = func(data, ifp, united)
        sections.append( section(title, comment, format_string, values) )
        ordering.append( order )

    io.write("%VERSION  VERSION_STAMP = V0001.000\n")#  DATE = 08/31/16 16:40:13\n")

    for order, sec in sorted(zip(ordering, sections)):
        io.write(sec)

def rst7(data, united, io):
    title = data.var["rnme"]
    io.write(title+'\n')
    atoms = helper.get_atoms(data, united)
    io.write(FortranWriter("I5,5E15.7").write([len(atoms),0]) + '\n')
    coords_nm = []
    [ coords_nm.extend(helper.get_coord(atom)) for atom in atoms ]
    coords = [ co * helper.NANOMETRE for co in coords_nm ]
    if DEBUG_VS_PDB:
        coords = [ round(co, 3) for co in coords ]

    numatoms = len(atoms)
    for i in range(numatoms):
        for j in range(i+1,numatoms):
            xi = helper.get_coord(atoms[i])
            xj = helper.get_coord(atoms[j])
            r = sqrt(sum((ii-jj)**2 for ii,jj in zip(xi,xj)))*10
            print("{} {} {}".format(i+1,j+1,r))

    io.write(FortranWriter("6F12.7").write(coords))
    io.write('\n')

def section(title, comment, format_string, values):
    header = section_header(title, comment, format_string)
    body = section_body(values, format_string)
    return header+body

def section_header(title, comment, format_string):
    flag = "%FLAG {}\n".format(title)
    nocomment = comment == amber_sections.NOCOMMENT
    nocomment = True
    com = "" if nocomment else "%COMMENT {}\n".format(title)
    fmt = "%FORMAT({})\n".format(format_string)
    return flag + com + fmt

def section_body(values, format_string):
    return FortranWriter(format_string).write(values) + '\n'

if __name__ == "__main__":
    dummy = dict()
    dummy['title'] = "hello world " * 100

    parm7(dummy, sys.stdout)

