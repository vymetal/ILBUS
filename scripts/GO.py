#!/usr/bin/python
import os, sys


os.system("./prepare_frames.py ndx frames_list")
dr = [i.strip() for i in open("DIRS", "r").readlines()]
f = open("done", "a")
f.close()

wd = os.getcwd()
for i in dr:
    fc = float(i)
    os.system("./prepare.py %f" % fc)

    try:
        f = open("done", "r")
        prev = f.readlines()[-1].strip()
        f.close()
    except:
        prev = ""

    os.system(
        "cd %s; pwd; for d in */; do cp $d/state.cpt ../$d/; done " % (wd + "/" + prev)
    )

    os.system("for i in `seq 1`; do . ./md_round ; done")
    os.system("./check_iteration.py")
    x = [i.strip() for i in open("noneq", "r").readlines() if i]
    cyc = 0
    while len(x) > 0:
        os.system("bash prolong")
        os.system("./check_iteration.py")
        x = [i.strip() for i in open("noneq", "r").readlines() if i]
        cyc += 1
        if cyc > 10:
            break

    if len(x) > 0:
        for iii in x:
            g = open(iii + "/no_conv", "w")
            g.close()

    os.mkdir(str(i))
    os.system("mv 0*/ " + i + "/")
    f = open("done", "a")
    f.write(str(i) + "\n")
    f.close()
