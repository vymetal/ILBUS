#!/usr/bin/python
import os, sys

fns = [i.strip() for i in open("frames_list", "r").readlines()]
x = []
try:
    x = [i.strip() for i in open("low_overlaps", "r").readlines()]
except:
    pass

fc = float(sys.argv[1])

wd = os.getcwd()
for ii, f in enumerate(fns):
    dr = "%06d" % ii
    os.system("mkdir %s" % dr)
    os.chdir(dr)
    if (not x) or (dr in x):
        os.system("cp -a ../amber14sb.ff ./")
        os.system("cp ../topol.top ./")
        os.system("cp %s ./" % f)
        os.system(
            "~/gromacs-2016/build/bin/gmx grompp -f ../run.mdp -c *.gro -maxwarn 2"
        )
        os.system("rm -r amber14sb.ff")
        os.system("rm topol.top")
        os.system("../gro_to_plumedfiles ../ndx ../frames_list %i %f" % (ii, fc))
    else:  # just create fake COLVAR
        f = open("COLVAR", "w")
        f.close()

    os.chdir(wd)
