#!/usr/bin/python
# join sets with the same landmarks but different force contants
# it writes only equilibrated data

import numpy
import scipy.stats

####EDIT THESE VARIBLES####
# force constants
drs = (
    "100000",
    "75000",
    "56250",
    "42188",
    "31641",
    "23731",
    "17798",
    "13348",
    "10011",
    "7508",
    "5631",
)
# number of landmarks
N = 29
###########################

EQ = 20  # equilibration times (in frames/samples)


def testks23(rms):
    n = len(rms)
    data21 = rms[: n // 2]
    data22 = rms[n // 2 :]
    data31 = rms[: n // 3]
    data32 = rms[n // 3 : 2 * n // 3]
    data33 = rms[2 * n // 3 :]
    lim = 0.1
    return (
        scipy.stats.ks_2samp(data21, data22)[1] > lim
        and scipy.stats.ks_2samp(data31, data32)[1] > lim
        and scipy.stats.ks_2samp(data31, data33)[1] > lim
        and scipy.stats.ks_2samp(data32, data33)[1] > lim
    )


def L(fn, my):
    l = open(fn, "r").readlines()
    l = [i.strip().split() for i in l if i.strip()[0] != "#"]
    # remove duplicate lines
    for j in range(len(l) - 1, 0, -1):
        if l[j] == l[j - 1]:
            del l[j]

    rs = l[EQ:]

    if len(rs) == 0:
        return rs

    rms = [float(r[2 + my]) for r in rs]

    for shift in range(1, len(rms)):
        cor = numpy.corrcoef(rms[:-shift], rms[shift:])[0, 1]
        if cor < 0.05:
            break
        if len(rms[::shift]) < 100:
            break
    rms = rms[::shift]

    def find_len(data):
        for cyc in range(len(data) - 2):
            s = testks23(data)
            if s:
                return len(data)
            data = data[1:]
        return 0

    ln = find_len(rms)
    start = (len(rms) - ln) * shift
    start = 0

    rs = rs[start:]
    return rs


fns = []
mys = []
for d in drs:
    for i in range(N):
        fns.append(d + "/" + "%06d" % i + "/COLVAR")
        mys.append(i)

import os, sys

try:
    os.mkdir("joined")
except:
    pass
for ii, cv in enumerate(fns):
    dr = "joined/%06d" % ii
    try:
        os.mkdir(dr)
    except:
        pass
    o = open(dr + "/COLVAR", "w")
    f = L(cv, mys[ii])
    print(cv, len(f))
    for line in f:
        out = line[:]
        for k in range(len(drs) - 1):
            out.extend(line[2:])
        o.write(" ".join(out) + "\n")
    o.close()

    # copy other files#
    src = cv.replace("COLVAR", "rama.dat")
    dst = dr + "/"
    os.system("cp %s %s" % (src, dst))
