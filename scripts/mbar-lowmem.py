#!/usr/bin/python
import sys, numpy, math, pylab as plt, random

EQUILIB0 = 0
RT = 8.314e-3 * 300
CUTOFF = 1e-6  # convergence cutoff
savedir = "./"
fake_times = False


def load(fn, my):
    data = open(fn).readlines()

    ndata = []
    k = 0
    while k < len(data):
        if "#" in data[k]:
            if data[k].find("#!") > -1:
                k += 1  # keep the value after exchange, because it matches the concatenated trajectory
            elif data[k].find("#EXCHANGE") > -1:
                del ndata[-1]  # get rid of the last frame before exchange
                k += EQUILIB + 1
            else:
                k += 1
        else:
            ndata.append(data[k])
            k += 1
    data = [i.strip().split() for i in ndata]
    tms = {}
    ndata = []

    if fake_times:
        for ii in range(len(data)):
            data[ii][0] = str(ii)

    for line in data:
        if line[0] not in tms:
            ndata.append(line)
            tms[line[0]] = True

    data = ndata

    x = [[float(j) for j in i[2:]] for i in data]
    bias = [float(i[1]) for i in data]
    tm = [i[0] for i in data]  # time stamp

    if len(x) > 0:
        A = numpy.array(
            [
                [
                    i[my] ** 2,
                ]
                for i in x
            ]
        )
        prm = numpy.linalg.lstsq(A, bias)[0]
        k = 2 * prm[0]
        x0 = 0  # -prm[1]/k
        print(k, x0)
        return x, k, x0, tm, bias
    else:  # it is empty - placeholder
        return None, None, None, None, None


sets = []
ks = []
x0s = []

maxframes = 100000000
excluded = ""
try:
    maxframes = int(sys.argv[1])
    excluded = sys.argv[1]
except:
    pass
print("maxframes=", maxframes)

idfs = []
bias = []
sims = []
mask = []
for i in sys.argv[1:]:
    if i not in (excluded,):

        sim = i.split("/")[-2]
        indx = int(sim)  # number of simulation

        s, k, x0, tm, bi = load(i, indx)
        if s != None:
            mask.append(True)
            sims.append(sim)
            s = s[
                EQUILIB0 : EQUILIB0 + maxframes
            ]  # discard initial frames for equilibration
            tm = tm[EQUILIB0 : EQUILIB0 + maxframes]
            bi = bi[EQUILIB0 : EQUILIB0 + maxframes]
            print("len(s)=", len(s))
            ls = len(s)
            ks.append(k)
            x0s.append(x0)
            idfs.extend([sim + " " + t for t in tm])
            bias.extend(bi)

            s = numpy.array(s)
            numpy.save(savedir + "set_" + sim + ".npy", s)
            print(s.shape)
        else:
            print("Zero length, neglecting")
            mask.append(False)
lens = []

for sim in sims:
    s = numpy.load(savedir + "set_" + sim + ".npy")
    lens.append(s.shape[0])
    s = s[:, mask]  # remove extra/empty columns
    v = 0.5 * numpy.multiply(numpy.array(ks)[None, :], numpy.square(s - x0s))
    expv = numpy.exp(-v / RT)
    numpy.save(savedir + "set_" + sim + ".npy", expv)  # replace the original data
    del v
    del expv

nB = len(sims)
Z = numpy.ones(nB)
iZ = numpy.reciprocal(Z)
for cyc in range(5000):
    Zold = Z.copy()
    Z = numpy.zeros(nB)
    w = []
    for sim in sims:
        expv = numpy.load(savedir + "set_" + sim + ".npy")
        ww = numpy.reciprocal(
            numpy.sum(numpy.multiply(numpy.multiply(lens, iZ)[None, :], expv), axis=1)
        )
        Z += numpy.sum(numpy.multiply(expv, ww[:, None]), axis=0)
        w.extend(ww)
        del expv
    w = numpy.array(w)
    print("\t\t\t", sum(w))
    Z = Z / numpy.sum(Z)
    w = w / numpy.sum(w)
    iZ = numpy.reciprocal(Z)
    eps = numpy.sum(numpy.square(numpy.log(numpy.divide(Z, Zold))), dtype=numpy.double)
    print(eps, max(w) / (1.0 / w.shape[0]))
    if eps < CUTOFF:
        break

print("...Effective sample sizes...")
k = 0
for s, l in zip(sims, lens):
    ww = w[k : k + l]
    neff = numpy.sum(ww) ** 2 / numpy.sum(numpy.square(ww))
    k += l
    print(s, ":", neff)
print("sum Z=", numpy.sum(Z))
print("sum w=", numpy.sum(w))

out = open("uwham_ks.dat", "w")
for i in ks:
    out.write(str(i) + "\n")
out.close()

out = open("uwham_ns.dat", "w")
for i in lens:
    out.write(str(i) + "\n")
out.close()

out = open("uwham_iZ.dat", "w")
for i in iZ:
    out.write(str(i) + "\n")
out.close()

out = open("uwham_mask.dat", "w")
for i in mask:
    out.write(str(i) + "\n")
out.close()


o = 0
bias = numpy.array(bias)

fe = -RT * numpy.log(Z)
fe = fe - min(fe)
for i in fe:
    print(i)

suf = ""

out = open("fes" + suf + ".dat", "w")
for i in fe:
    out.write(str(i) + "\n")
out.close()
#
out = open("W" + suf + ".dat", "w")
for t, wh in zip(idfs, w):
    out.write(t + "\t" + str(wh) + "\n")
out.close()
