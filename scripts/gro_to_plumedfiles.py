#!/usr/bin/python
import sys

print(
    sys.argv[0],
    " [index file] [list of reference structures] [ my index ] [ force constant ]",
)

ndxfile = sys.argv[1]
myind = int(sys.argv[3])
fc = float(sys.argv[4])

# load atomic indices for rmsd calculation, starting from 1
ndx = []
for i in open(ndxfile, "r").readlines():
    i = i.strip().split()
    for x in i:
        ndx.append(int(x))

list0 = [i.strip() for i in open(sys.argv[2], "r").readlines()]
list0 = [i for i in list0 if i]

nms = []
for fn in list0:
    base = fn.split("/")[-1]
    base = base.split(".")[0]
    nm = base + ".pdb"
    nms.append(nm)

out = open("plumed.dat", "w")
out.write("UNITS LENGTH=nm TIME=fs\n")
out.write("WHOLEMOLECULES ENTITY0=1-62\n")  # edit these number for different molecule
for ii, nm in enumerate(nms):
    out.write("r%i: RMSD REFERENCE=%s TYPE=OPTIMAL\n" % (ii, "../frames/" + nm))
out.write("restraint: RESTRAINT ARG=r%i  KAPPA=%f  AT=0.0\n" % (myind, fc))
out.write(
    "PRINT ARG=restraint.bias,%s STRIDE=500 FILE=COLVAR\n"
    % (",".join(["r%i" % i for i in range(len(nms))]))
)
out.close()
