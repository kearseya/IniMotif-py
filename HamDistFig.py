from KmerKounter import identifier

from KmerKounter import kmercount
from KmerKounter import numofruns
from KmerKounter import hash2kmer
from KmerKounter import kmer2hash
from KmerKounter import revcompwanted
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import inputlist
from KmerKounter import startround
from KmerKounter import totaldict

from KmerKounter import top6all

from WebLogoMod import nmotifs
from WebLogoMod import removelist

import matplotlib.pyplot as plt
import numpy as np
import math
import random
#import operator

from adjustText import adjust_text

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
def hammer():
    for _ in range(numofruns + 1):
        khams.append({})
    for run in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            khams[run][k] = {}
            for n in range(1, nmotifs+1):
                khams[run][k][n] = {k: {} for k in range(k+1)}
                #print("run "+str(run)+" k "+str(k))
                #print(khams)
                if n == 1:
                    hconsensus = max(kmercount[run][k], key=lambda key: kmercount[run][k][key])
                else:
                    temptop6 = top6all[run][k].copy()
                    for j in temptop6:
                        if j in removelist[run][k][n]:
                            temptop6.remove(j)
                            if kmer2hash(revComp(hash2kmer(j, k))) in temptop6:
                                temptop6.remove(kmer2hash(revComp(hash2kmer(j, k))))
                    hconsensus = max(temptop6, key=lambda key: kmercount[run][k][key])
                #hconsensus = max(kmercount[run][k], key=lambda key: kmercount[run][k][key])
                consensus = hash2kmer(hconsensus, k)
                #print("con: "+str(consensus))
                for x in list(kmercount[run][k].keys()):
                    if x not in removelist[run][k][n]:
                        values = hash2kmer(x,k)
                        rvalues = revComp(values)
                        ham = hamming_distance(consensus, values)
                        rham = hamming_distance(consensus, rvalues)
                        if ham <= rham:
                            khams[run][k][n][ham].update({x:kmercount[run][k][x]})
                        if ham > rham:
                            khams[run][k][n][rham].update({x:kmercount[run][k][x]})
    return khams

hammer()
"""
def multihammer():
    for k in range (mink, maxk+1):
        hammer(k)
"""
#multihammer()
#print("khams")
#print(khams)

removedkmers = []

def findremovedkmers():
    global removedkmers
    for _ in range(0, numofruns+1):
        removedkmers.append({})
    for r in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            removedkmers[r][k] = {}
            for n in range(2, nmotifs+1):
                removedkmers[r][k][n] = 0
                for kmers in removelist[r][k][n]:
                    try:
                        removedkmers[r][k][n] += kmercount[r][k][kmers]
                    except:
                        continue

findremovedkmers()

#print("numremovedkmers")
#print(removedkmers)

#print("removelist")
#print(removelist)

"""
top12t = []
top12 = []

def top12maker(numofruns):
    for _ in range(0, numofruns+1):
        top12t.append({})
        top12.append({})
    for x in range(1, numofruns+1):
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for k in range(mink, maxk+1):
            top12t[x][k] = list(sorted(kmercount[x][k].items(), key=lambda x: x[1], reverse=True))[:12]
            top12[x][k] = [j[0] for j in top12t[x][k]]

top12maker(numofruns)
#print(kmercount)
#print(top12)

top6all = []
"""

"""
def top6(run, k):
    top6s = []
    #keys = top12[run][k]
    topvalpos = 0
    while len(top6s) <= 11:
        next = top12[run][k][topvalpos]
        nkmer = hash2kmer(next, k)
        nrkmer = revComp(nkmer)
        nhrkmer = kmer2hash(nrkmer)
        if next not in top6s:
            top6s.append(next)
            #print("k "+str(k)+" len "+str(len(top6s))+" kmer "+str(next))
        if next == nhrkmer:
            if len(top6s) <= 11:
                top6s.append(nhrkmer)
        try:
            if kmercount[run][k][nhrkmer] > 0:
                if nhrkmer not in top6s:
                    if len(top6s) <= 11:
                        top6s.append(nhrkmer)
                        #print("k "+str(k)+" len "+str(len(top6s))+" kmer "+str(nhrkmer))
        except:
            if len(top6s) <= 11:
                top6s.append(next)
        topvalpos += 1
    return top6s
"""

"""
def top6(run, k):
    top6us = []
    top6ts = []
    top6ds = {}
    top6s = []
    topvalpos = 0
    while len(top6us) <= 11:
        next = top12[run][k][topvalpos]
        nkmer = hash2kmer(next, k)
        nrkmer = revComp(nkmer)
        nhrkmer = kmer2hash(nrkmer)
        if next not in top6us:
            top6us.append(next)
            #print("k "+str(k)+" len "+str(len(top6s))+" kmer "+str(next))
        if next == nhrkmer:
            if len(top6us) <= 11:
                top6us.append(nhrkmer)
        try:
            if kmercount[run][k][nhrkmer] > 0:
                if nhrkmer not in top6us:
                    if len(top6us) <= 11:
                        top6us.append(nhrkmer)
                        #print("k "+str(k)+" len "+str(len(top6s))+" kmer "+str(nhrkmer))
        except:
            if len(top6us) <= 11:
                top6us.append(next)
        topvalpos += 1

    for i in range(len(top6us)):
        if i//2 == (i+1)//2:
            top6ts.append((top6us[i], top6us[i+1]))

    for x in top6ts:
        try:
            first = kmercount[run][k][x[0]]
        except:
            first = 0
        if x[0] != x[1]:
            try:
                second = kmercount[run][k][x[1]]
            except:
                second = 0
        if x[0] == x[1]:
            second = 0
        top6ds[x] = first+second

    sort = sorted(top6ds.items(), key=operator.itemgetter(1), reverse=True)

    for z in range(len(sort)):
        top6s.append(sort[z][0][0])
        top6s.append(sort[z][0][1])

    return top6s




def top6maker(numofruns):
    for _ in range(0, numofruns+1):
        top6all.append({})
    for x in range(1, numofruns+1):
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for k in range(mink, maxk+1):
            top6all[x][k] = top6(x, k)

top6maker(numofruns)
#print("top6all")
#print(top6all)
"""


def top6plotx(run, p, k, n):
    x = []
    keys = [*khams[run][k][n][p].keys()]
    vals = []
    if n == 1:
        top6s = top6all[run][k].copy()
    else:
        top6s = top6all[run][k].copy()
        for j in top6s:
            if j in removelist[run][k][n]:
                top6s.remove(j)
                rj = kmer2hash(revComp(hash2kmer(j, k)))
                if rj in top6s:
                    top6s.remove(rj)
    #print("TOP6Sx")
    #print(top6s)
    for i in top6s:
        if i in keys:
            vals.append(i)
    for _ in range(len(vals)):
        nval = float(p) + random.uniform(-0.3, 0.3)
        x.append(nval)
    return x

def top6ploty(run, x, k, n):
    y = []
    values = [*khams[run][k][n][x].keys()]
    vals = []
    if n == 1:
        top6s = top6all[run][k].copy()
    else:
        top6s = top6all[run][k].copy()
        for j in top6s:
            if j in removelist[run][k][n]:
                top6s.remove(j)
                rj = kmer2hash(revComp(hash2kmer(j, k)))
                if rj in top6s:
                    top6s.remove(rj)
    #print("TOP6Sy")
    #print(top6s)
    for i in top6s:
        if i in values:
            vals.append(i)
    for i in vals:
        val = khams[run][k][n][x][i]
        nval = float(val) + random.uniform(-0.1, 0.1)
        y.append(nval)
    return y



def top6split(run, x, k, n):
    values = [*khams[run][k][n][x].keys()]
    vals = []
    if n == 1:
        top6s = top6all[run][k].copy()
    else:
        top6s = top6all[run][k].copy()
        for j in top6s:
            if j in removelist[run][k][n]:
                top6s.remove(j)
                rj = kmer2hash(revComp(hash2kmer(j, k)))
                if rj in top6s:
                    top6s.remove(rj)

    for i in top6s[::2]:
        if i in values:
            vals.append(i)
            rk = kmer2hash(revComp(hash2kmer(i,k)))
            vals.append(rk)
            try:
                if rk not in khams[run][k][n][x]:
                    khams[run][k][n][x][rk] = 0
            except:
                continue
    return vals


def top6splitter(run, k, n):
    split = []
    for _ in range(numofruns + 1):
        split.append({})
    for x in khams[run][k][n]:
        top6sp = top6split(run, x, k, n)
        split[run][x] = top6sp
    return split



def xaxismaker(run, p, k, n):
    x = []
    keys = [*khams[run][k][n][p].keys()]
    if n == 1:
        top6s = top6all[run][k].copy()
    else:
        top6s = top6all[run][k].copy()
        for j in top6s:
            if j in removelist[run][k][n]:
                top6s.remove(j)
                rj = kmer2hash(revComp(hash2kmer(j, k)))
                if rj in top6s:
                    top6s.remove(rj)
    for i in top6s:
        if i in keys:
            keys.remove(i)
    size = len(keys)
    for _ in range(size):
        nval = float(p) + random.uniform(-0.3, 0.3)
        x.append(nval)
    return x

def yaxismaker(run, x, k, n):
    y = []
    values = [*khams[run][k][n][x].keys()]
    if n == 1:
        top6s = top6all[run][k].copy()
    else:
        top6s = top6all[run][k].copy()
        for j in top6s:
            if j in removelist[run][k][n]:
                top6s.remove(j)
                rj = kmer2hash(revComp(hash2kmer(j, k)))
                if rj in top6s:
                    top6s.remove(rj)
    for i in top6s:
        if i in values:
            values.remove(i)
    for i in values:
        val = khams[run][k][n][x][i]
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
labelcolours = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

#print("khams")
#print(khams)

def scatter(run, k, n):
    #top6(run, k)
    split = top6splitter(run, k, n)
    #print("split")
    #print(split)
    labels = []
    handels = []
    #labelcolours = []
    #greys = [0.6, 0.75, 0.9]
    c = 0
    cc = 0
    texts = []

    for x in khams[run][k][n]:

        #print("run: "+str(run)+", k: "+str(k)+", n: "+str(n))
        #print("X: "+str(x))

        xaxis = xaxismaker(run, x, k, n)
        yaxis = yaxismaker(run, x, k, n)

        top6x = top6plotx(run, x, k, n)
        #print("top6x")
        #print(top6x)
        top6y = top6ploty(run, x, k, n)
        #print("top6y")
        #print(top6y)

        plt.scatter(xaxis, yaxis, color = '0.75', alpha=0.7, s=1)

        count = len(top6x)

        for i in range(0, count):

            plt.scatter(top6x[i], top6y[i], label=split[run][x][i], color=colours[c], s=3)
            labels.append((str(hash2kmer(split[run][x][i],k)+' '+str(khams[run][k][n][x][(split[run][x][i])]))))
            c += 1



        for p, j in enumerate(split[run][x]):
            #print("split[run][x]")
            #print(split[run][x])
            #print("p: "+str(p))
            #print("j: "+str(j))
            #print("top6x: "+str(top6x[p]))
            #print("top6y: "+str(top6y[p]))

            xp = top6x[p]
            yp = top6y[p]

            texts.append(plt.text(float(xp), float(yp-0.3), str(hash2kmer(j,k) + ' ' + str(khams[run][k][n][x][j])),color=colours[cc], fontsize=5))

            #if (p % 2)+x == x:
                #plt.annotate(str(hash2kmer(j,k) + ' ' + str(khams[run][k][n][x][j])), (float(xp), float(yp)+0.6), color=colours[cc], fontsize=5)
            #if (p % 2)+x != x:
                #plt.annotate(str(hash2kmer(j,k) + ' ' + str(khams[run][k][n][x][j])), (float(xp), float(yp)-0.6), color=colours[cc], fontsize=5)

            cc += 1

    adjust_text(texts, lw=0.5)
    #for z in colours[::2]:
        #labelcolours.append(z)

    leg = plt.legend(labels[::2], fontsize=7)
    if n == 1:
        top6s = top6all[run][k].copy()
    else:
        top6s = top6all[run][k].copy()
        for j in top6s:
            if j in removelist[run][k][n]:
                top6s.remove(j)
                rj = kmer2hash(revComp(hash2kmer(j, k)))
                if rj in top6s:
                    top6s.remove(rj)

    for i in range(0,(len(top6s)//2)):

        leg.legendHandles[i].set_color(labelcolours[i])
        leg.legendHandles[i]._sizes = [8]

    plt.xticks(np.arange(0, len(khams[run][k][n]), step=1))

    plt.xlabel("Hamming distance")
    plt.ylabel("Kmer count")
    if n == 1:
        plt.title("Run number:"+' '+str(run+(startround-1))+'\n'+'K: '+str(k)+'\n'+'Total kmers: '+str(totaldict[run][k]))
    else:
        plt.title("Run number:"+' '+str(run+(startround-1))+'\n'+'K: '+str(k)+'\n'+'Total kmers: '+str(totaldict[run][k]-removedkmers[run][k][n]))
    #plt.show()
    plt.savefig("figures/"+str(identifier)+"/hamming_distance/hamdist_"+str(identifier)+"_"+str(run+(startround-1))+"_"+str(k)+"_"+str(n), dpi=600)
    plt.close()

def plotrange(runs):
    for r in range(1, runs+1):
        #mink = int(inputlist[((r-1)*5)+7])
        #maxk = int(inputlist[((r-1)*5)+8])
        for k in range(mink, maxk+1):
            for n in range(1, nmotifs+1):
                scatter(r, k, n)

plotrange(numofruns)
