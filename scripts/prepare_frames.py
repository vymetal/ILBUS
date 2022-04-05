#!/usr/bin/python
import sys, os

print(sys.argv[0], " [index file] [list of reference structures]")

ndxfile = sys.argv[1]

# load atomic indices for rmsd calculation, starting from 1
ndx = []
for i in open(ndxfile, "r").readlines():
    i = i.strip().split()
    for x in i:
        ndx.append(int(x))

list0 = [i.strip() for i in open(sys.argv[2], "r").readlines()]
list0 = [i for i in list0 if i]


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
            data.append((resnum, resname, aname, anum, xyz))
    return data


nms = []
os.system("mkdir frames")
for fn in list0:
    base = fn.split("/")[-1]
    base = base.split(".")[0]
    atoms = load_gro(fn, ndx)
    nm = base + ".pdb"
    nms.append(nm)

    out = open("frames/" + nm, "w")
    for a in atoms:
        resnum, resname, aname, anum, xyz = a
        out.write(
            "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
            % (
                "ATOM",
                anum,
                aname,
                resname,
                " ",
                resnum,
                " ",
                xyz[0] * 10.0,
                xyz[1] * 10.0,
                xyz[2] * 10.0,
                1.0,
                1.0,
                "",
            )
        )
    out.close()
