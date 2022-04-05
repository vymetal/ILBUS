#!/usr/bin/python
import os, sys


dr = [i.strip() for i in open("DIRS", "r").readlines()]
f = open("done", "w")
f.close()

wd = os.getcwd()
for i in dr:
    fc = float(i)
    os.system("mv %s/* ./" % i)
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

    os.system("mv 0*/ " + i + "/")
    f = open("done", "a")
    f.write(str(i) + "\n")
    f.close()
