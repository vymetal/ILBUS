#!/usr/bin/python

import sys, os, pylab as plt, numpy
import scipy.stats

###EDIT THIS VARIABLE############
NN = 30  # number of landmarks
###OTHER PARAMETERS OF ILBUS ####
FR = 100000000  # maximal number frames to analyze
EQ = 20  # number of initial frames to discard as equilibration
OL = 0.33  # the requested (total) overlap between a window and the remaining windows
MIN_LEN = 500  # minimal length of the simulation (in frames)
#################################


dirs = open("done", "r").readlines()
dirs = [i.strip() for i in dirs]
if len(sys.argv) > 1:
    dirs.append(sys.argv[1])
else:
    dirs.append("./")

dirs = dirs[::-1]


def testks23(rms):
    n = len(rms)
    data21 = rms[: n // 2]
    data22 = rms[n // 2 :]
    data31 = rms[: n // 3]
    data32 = rms[n // 3 : 2 * n // 3]
    data33 = rms[2 * n // 3 :]
    lim = 0.1  # 0.05
    return (
        scipy.stats.ks_2samp(data21, data22)[1] > lim
        and scipy.stats.ks_2samp(data31, data32)[1] > lim
        and scipy.stats.ks_2samp(data31, data33)[1] > lim
        and scipy.stats.ks_2samp(data32, data33)[1] > lim
    )


def L3(fn, my):
    l = open(fn, "r").readlines()
    l = [i.strip().split() for i in l if i.strip()[0] != "#"]
    # remove duplicate lines
    for j in range(len(l) - 1, 0, -1):
        if l[j] == l[j - 1]:
            del l[j]

    rs = [[float(j) for j in i[2:]] for i in l[EQ : FR + EQ]]
    if len(rs) == 0:
        return None

    rms = [r[my] for r in rs]

    for shift in range(1, len(rms)):
        cor = numpy.corrcoef(rms[:-shift], rms[shift:])[0, 1]
        if cor < 0.05:
            break
        if len(rms[::shift]) < 100:
            break
    #   print('l3:',shift,cor)
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

    rs = rs[start:]
    nfr = len(rs)

    n = NN  # len(rs[0])
    mat = numpy.zeros(n)
    for r in rs:
        for ii, i in enumerate(r):
            if i <= r[my]:
                mat[ii] += 1.0
    mat = mat / len(rs)
    lim = 0.1  # 0.025
    # export overlaps
    f = open("connect", "a")
    v = []
    for i in range(NN):
        if mat[i] > lim and i != my:
            v.append([mat[i], i])
    if v:
        v.sort(reverse=True)
    f.write("%d->" % my)
    for vv in v:
        f.write(" %d(%f)" % (vv[1], vv[0]))
    f.write("\n")
    f.close()
    # lim=len(rs)*0.025  #5%
    mat[mat < lim] = 0
    mat[mat > 0] = 1.0
    return mat


kk = 0
dr = dirs[0]
d = []
means = []
maxs = []
mat = []
kstest = []
mat = numpy.zeros([NN, NN])
for dr in dirs:
    for k in range(NN):
        m = L3(dr + "/%06d/COLVAR" % (k), k)
        if type(m) != type(None):
            mat[k, :] += m

print(mat)

mat = mat + mat.T

n = mat.shape[0]
vis = [False for i in range(n)]
tvis = [False for i in vis]
xcyc = 0
groups = []
while False in tvis:
    vis = [False for i in range(n)]
    ii = tvis.index(False)
    vis[ii] = True
    for cyc in range(2 * n):
        for i in range(n):
            if vis[i]:
                for j in range(n):
                    if mat[i, j] > 0:
                        vis[j] = True
        if vis.count(True) == n:
            break
    tvis = [a or b for a, b in zip(vis, tvis)]
    xcyc += 1
    print("group:", vis.count(True), [ii for ii, i in enumerate(vis) if i])
    groups.append([ii for ii, i in enumerate(vis) if i])
if xcyc == 1:
    print(True)
else:
    print(False)

rm = numpy.loadtxt("rmsd_mat")

### find candidates for extra landmarks
if len(groups) > 1:
    for i1 in range(len(groups)):
        for i2 in range(i1 + 1, len(groups)):
            g1 = groups[i1]
            g2 = groups[i2]
            minr = [g1[0], g2[0]]
            for j in g1:
                for k in g2:
                    if rm[j, k] < rm[minr[0], minr[1]]:
                        minr = [j, k]
            print(i1, i2, ":", minr[0], minr[1], rm[minr[0], minr[1]])
