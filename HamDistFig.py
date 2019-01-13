from KmerKounter import kmercount
from KmerKounter import numofruns
from KmerKounter import hash2kmer
from KmerKounter import kmer2hash
from KmerKounter import revcompwanted
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import runlists

import matplotlib.pyplot as plt
import numpy as np
import math
import random



revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



khams = []
def hammer(k):
    for _ in range(numofruns + 1):
        khams.append({})
    for j in range(1,len(kmercount)):
        khams[j][k] = {k: {} for k in range(k+1)}
        hconsensus = max(kmercount[j][k], key=lambda key: kmercount[j][k][key])
        consensus = hash2kmer(hconsensus, k)
        for x in list(kmercount[j][k].keys()):
            values = hash2kmer(x,k)
            rvalues = revComp(values)
            ham = hamming_distance(consensus, values)
            rham = hamming_distance(consensus, rvalues)
            if ham <= rham:
                khams[j][k][ham].update({x:kmercount[j][k][x]})
            if ham > rham:
                khams[j][k][rham].update({x:kmercount[j][k][x]})
    return khams

def multihammer(mink, maxk):
    for k in range (mink, maxk+1):
        hammer(k)

multihammer(mink, maxk)
print(khams)


def top6(run, k):
    top6s = []
    keys = list(kmercount[run][k].keys())
    topvalpos = 0
    while len(top6s) <= 11:
        next = keys[topvalpos]
        nkmer = hash2kmer(next, k)
        nrkmer = revComp(nkmer)
        nhrkmer = kmer2hash(nrkmer)
        if next not in top6s:
            top6s.append(next)
        if nhrkmer not in top6s:
            top6s.append(nhrkmer)

        topvalpos += 1
    return top6s

def top6plotx(run, p, k):
    x = []
    keys = [*khams[run][k][p].keys()]
    vals = []
    top6s = top6(run, k)
    for i in top6s:
        if i in keys:
            vals.append(i)
    for _ in range(len(vals)):
        nval = float(p) + random.uniform(-0.3, 0.3)
        x.append(nval)
    return x

def top6ploty(run, x, k):
    y = []
    values = [*khams[run][k][x].keys()]
    vals = []
    top6s = top6(run, k)
    for i in top6s:
        if i in values:
            vals.append(i)
    for i in vals:
        val = khams[run][k][x][i]
        nval = float(val) + random.uniform(-0.1, 0.1)
        y.append(nval)
    return y



def top6split(run, x, k):
    values = [*khams[run][k][x].keys()]
    vals = []
    top6s = top6(run, k)
    for i in top6s:
        if i in values:
            vals.append(i)
    return vals


def top6splitter(run, k):
    split = []
    for _ in range(numofruns + 1):
        split.append({})
    for x in khams[run][k]:
        top6sp = top6split(run, x, k)
        split[run][x] = top6sp
    return split



def xaxismaker(run, p, k):
    x = []
    keys = [*khams[run][k][p].keys()]
    top6s = top6(run, k)
    for i in top6s:
        if i in keys:
            keys.remove(i)
    size = len(keys)
    for _ in range(size):
        nval = float(p) + random.uniform(-0.3, 0.3)
        x.append(nval)
    return x

def yaxismaker(run, x, k):
    y = []
    values = [*khams[run][k][x].keys()]
    top6s = top6(run, k)
    for i in top6s:
        if i in values:
            values.remove(i)
    for i in values:
        val = khams[run][k][x][i]
        nval = float(val) + random.uniform(-0.1, 0.1)
        y.append(nval)
    return y


"""
FOR TESTING
topx = []

def converter(run, k):
    for x in top6:
        topx.append(hash2kmer(x,6))

converter(1,6)
print(topx)
"""

colours = ['C0', 'C0', 'C1', 'C1', 'C2', 'C2', 'C3', 'C3', 'C4', 'C4', 'C5', 'C5', 'C6', 'C6', 'C7', 'C7']


def scatter(run, k):
    top6(run, k)
    split = top6splitter(run, k)
    labels = []
    handels = []
    labelcolours = []
    #greys = [0.6, 0.75, 0.9]
    c = 0
    cc = 0

    for x in khams[run][k]:

        xaxis = xaxismaker(run, x, k)
        yaxis = yaxismaker(run, x, k)

        top6x = top6plotx(run, x, k)
        top6y = top6ploty(run, x, k)

        plt.scatter(xaxis, yaxis, color = '0.75', alpha=0.7)

        count = len(top6x)

        for i in range(0, count):
            plt.scatter(top6x[i], top6y[i], label=split[run][x][i], color=colours[c])
            labels.append((str(hash2kmer(split[run][x][i],k)+' '+str(khams[run][k][x][(split[run][x][i])]))))
            c += 1

        for p, j in enumerate(split[run][x]):
            #print(p)
            xp = top6x[p]
            yp = top6y[p]
            #print(j)
            #print(xp)
            #print(yp)
            plt.annotate(str(hash2kmer(j,k) + ' ' + str(khams[run][k][x][j])), (float(xp), float(yp)), color=colours[cc])
            cc += 1

    for z in colours[::2]:
        labelcolours.append(z)

    print(handels)
    print(labels)
    leg = plt.legend(labels[::2])
    for i in range(0,6):
        leg.legendHandles[i].set_color(labelcolours[i])

    plt.xlabel("Hamming distance")
    plt.ylabel("Kmer count")
    plt.title("Run number:"+' '+str(run)+'\n'+'K: '+str(k)+'\n'+'Total kmers: '+str(len(runlists[run][k])))
    plt.show()

def plotrange(runs, mink, maxk):
    for r in range(1,runs+1):
        for k in range(mink, maxk+1):
            scatter(r, k)

plotrange(numofruns, mink, maxk)
