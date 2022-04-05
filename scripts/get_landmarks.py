#!/usr/bin/python
import numpy, math, sys, hashlib

print(sys.argv[0], "[index file] [frame_list] [min RMSD]")

minrmsd = float(sys.argv[3])
ndxfile = sys.argv[1]

# load atomic indices for rmsd calculation, starting from 1
ndx = []
for i in open(ndxfile, "r").readlines():
    i = i.strip().split()
    for x in i:
        ndx.append(int(x))

list1 = [i.strip() for i in open(sys.argv[2], "r").readlines()]
list1 = [i for i in list1 if i]

print(len(ndx))


def load_gro(fn, filt=None):
    print("loading %s..." % fn)
    f = open(fn, "r")
    lines = f.readlines()
    f.close()
    data = []
    for line in lines[2:-1]:
        resnum = int(line[0:5])
        resname = line[5:10].strip()
        aname = line[10:15].strip()
        anum = int(line[15:20])
        xyz = (float(line[20:28]), float(line[28:36]), float(line[36:44]))
        if not (filt) or (anum in filt):
            data.append(xyz)
    return data


CALC = 0


def rmsd(ref, s):
    com1 = numpy.mean(ref, axis=0)
    com2 = numpy.mean(s, axis=0)
    sel1 = numpy.array(ref) - com1
    sel2 = numpy.array(s) - com2
    V, S, Wt = numpy.linalg.svd(numpy.dot(numpy.transpose(sel2), sel1))
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
    if reflect < -0.1:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    U = numpy.dot(V, Wt)
    new = numpy.dot(sel2, U)
    dx = (new - sel1).flatten()
    n = len(dx) / 3
    global CALC
    CALC += 1
    return math.sqrt(numpy.dot(dx, dx) / n)


#####################################
refs = []
frnames = []

for fn in list1:
    xx = load_gro(fn, ndx)
    if not (refs):
        refs.append(xx)
        frnames.append(fn)
    rmss = [rmsd(xx, x0) for x0 in refs]
    if min(rmss) > minrmsd:
        refs.append(xx)
        frnames.append(fn)

print(len(refs), " frames")
f = open("frames_list", "w")
for i in frnames:
    f.write(i + "\n")
f.close()
