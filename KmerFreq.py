from KmerKounter import identifier

from KmerKounter import kmercount
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import kmer2hash
from KmerKounter import hash2kmer
from KmerKounter import revnuc
from KmerKounter import revComp

from KmerKounter import numofruns
from HamDistFig import top6all

import matplotlib.pyplot as plt
from matplotlib.transforms import TransformedBbox, Bbox
import numpy as np
import math



def makexaxis():
    xaxis = ([])
    for i in range(numofruns):
        xaxis.append(i)
    return xaxis

"""
def top6(k):
    top6s = []
    keys = list(kmercount[(numofruns)][k].keys())
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
"""

#print(top6all)
"""
def colours1(k):
    colours = []
    top6s = top6all[numofruns][k]
    consensusnum = kmercount[numofruns][k][top6s[0]]
    for i in range(0, 4**k):
        if kmercount[numofruns][k][i] not in top6s:
            if kmercount[numofruns][k][i] >= consensusnum//5:
                colours.append(i)
    return colours
"""


def colours1():
    #colours = []
    global colours
    colours = {}
    for k in range(mink, maxk+1):
        top6s = top6all[numofruns][k]
        #consensusnum = kmercount[numofruns][k][top6s[0]]
        colours[k] = list(sorted(kmercount[numofruns][k].items(), key=lambda x: x[1], reverse=True))[:62]
        colours[k] = [j[0] for j in colours[k]]
        for i in top6s:
            if i in colours[k]:
                colours[k].remove(i)
    return colours

colours1()

"""
def makexaxis2():
    xaxis = ([])
    for i in range(-1, len(kmercount)+1):
        xaxis.append(i)
    return xaxis
"""

def makeyaxis1a(i, k):
    yaxis1a = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    #colours = colours[k]
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        if i not in top6s or colours[k]:
            num =  kmercount[l][k][i]
            rnum = kmercount[l][k][rkmer]
            f = (num+rnum)/total
            freq = (f/(1-f))
            yaxis1a.append(math.log10(freq))
    return yaxis1a



def makeyaxis1b(i, k):
    yaxis1b = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    #colours = colours[k]
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        if i in colours[k]:
            num =  kmercount[l][k][i]
            rnum = kmercount[l][k][rkmer]
            f = (num+rnum)/total
            freq = (f/(1-f))
            yaxis1b.append(math.log10(freq))
    return yaxis1b


def makeyaxis1c(i, k):
    yaxis1c = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        if i in top6s:
            num =  kmercount[l][k][i]
            rnum = kmercount[l][k][rkmer]
            f = (num+rnum)/total
            freq = (f/(1-f))
            yaxis1c.append(math.log10(freq))
    return yaxis1c



def makeyaxis2a(i, k):
    yaxis2a = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    #colours = colours[k]
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        if i not in top6s or colours[k]:
            num =  kmercount[l][k][i]
            rnum = kmercount[l][k][rkmer]
            yaxis2a.append((num+rnum)/total)
    return yaxis2a


def makeyaxis2b(i, k):
    yaxis2b = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    #colours = colours[k]
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        if i in colours[k]:
            num =  kmercount[l][k][i]
            rnum = kmercount[l][k][rkmer]
            yaxis2b.append((num+rnum)/total)
    return yaxis2b


def makeyaxis2c(i, k):
    yaxis2c = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        if i in top6s:
            num =  kmercount[l][k][i]
            rnum = kmercount[l][k][rkmer]
            yaxis2c.append((num+rnum)/total)
    return yaxis2c



def makeyaxis3(k):
    yaxis3 = ([])
    for l in range(1, (numofruns+1)):
        total = sum(kmercount[l][k].values())
        yaxis3.append(total)
    return yaxis3



def grapher(k):

    xaxis = makexaxis()
    last = (len(xaxis)-1)
    top6s = top6all[numofruns][k]
    #colourp = colours[k]

    fig = plt.figure(figsize=(10,10))
    grid = plt.GridSpec(2,3,wspace=0.4,hspace=0.3)

    top = fig.add_subplot(grid[:-1,:])
    top.set_xlabel("SELEX round")
    top.set_ylabel("log(f/(1-f))")
    top.set_title("Kmer frequency"+', K: '+str(k))
    top.set_xlim([-0.5, (numofruns+1)])
    top.set_xticks(np.linspace(-0.5, (numofruns+1), num=((2*(numofruns+1))+2), endpoint=True))

    bottom = fig.add_subplot(grid[-1,:-1])
    bottom.set_xlabel("SELEX round")
    bottom.set_ylabel("f = (kmer/total)")
    bottom.set_title("Kmer frequency")
    bottom.set_xlim([-0.5, numofruns])
    bottom.set_xticks(range(0,(numofruns+1)))

    bar = fig.add_subplot(grid[-1,-1:])
    bar.set_xlabel("SELEX round")
    bar.set_ylabel("Total kmers")
    bar.set_title("Kmer total distribution")
    bar.set_xticks(makexaxis())

    colourslist = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

    for l in range(1, (numofruns+1)):
        #total = sum(kmercount[l][k].values())
        for i in kmercount[l][k]:
            num = kmercount[l][k][i]
            rkmer = kmer2hash(revComp(hash2kmer(i,k)))
            try:
                rnum = kmercount[l][k][rkmer]
                if num > rnum:
                    if i in top6s:
                        yaxis1c = makeyaxis1c(i, k)
                        top.plot(xaxis, yaxis1c, color=colourslist[(top6s.index(i)//2)], linewidth=2, marker="s", markevery=None, zorder=int(kmercount[l][k][i]))
                    if i in colours[k]:
                        yaxis1b = makeyaxis1b(i, k)
                        top.plot(xaxis, yaxis1b, linewidth=1,  marker=".", markevery=None, zorder=int(kmercount[l][k][i]))
                    if i not in top6s or colours[k]:
                        yaxis1a = makeyaxis1a(i, k)
                        top.plot(xaxis, yaxis1a, color = '0.75', linestyle='--', linewidth=0.5, marker="x", markevery=None, alpha=0.5, zorder=0)
            except:
                continue


    for l in range(1, (numofruns+1)):
        #total = sum(kmercount[l][k].values())
        for i in kmercount[l][k]:
            num = kmercount[l][k][i]
            rkmer = kmer2hash(revComp(hash2kmer(i,k)))
            try:
                rnum = kmercount[l][k][rkmer]
                if num > rnum:
                    if i in top6s:
                        yaxis2c = makeyaxis2c(i, k)
                        bottom.plot(xaxis, yaxis2c, color=colourslist[(top6s.index(i)//2)], linewidth=2, marker="s", markevery=None, zorder=int(kmercount[l][k][i]))
                    if i in colours[k]:
                        yaxis2b = makeyaxis2b(i, k)
                        bottom.plot(xaxis, yaxis2b, linewidth=1,  marker=".", markevery=None, zorder=int(kmercount[l][k][i]))
                    if i not in top6s or colours[k]:
                        yaxis2a = makeyaxis2a(i, k)
                        bottom.plot(xaxis, yaxis2a, color = '0.75', linestyle='--', linewidth=0.5, marker="x", markevery=None, alpha=0.8, zorder=0)
            except:
                continue


    ymint, ymaxt = top.get_ylim()
    ypost = np.linspace(ymint, ymaxt, num=20, endpoint=True)

    yminb, ymaxb = bottom.get_ylim()
    yposb = np.linspace(yminb, ymaxb, num=20, endpoint=True)

    top6labels = []
    for n, i in enumerate(top6s[::2]):
        p = (n//2)
        yaxis1c = makeyaxis1c(i, k)
        yaxis2c = makeyaxis2c(i, k)
        top.annotate((str(n+1)), (xaxis[last], yaxis1c[last]), (xaxis[last]+0.1, ypost[-(n+2)]), size=10, fontname='monospace', weight='bold', arrowprops=dict(color=colourslist[n], shrink=0.05, width=0.05, headwidth=0.4), color=colourslist[n])
        top6labels.append(str(n+1)+". "+str(hash2kmer(i,k))+" / "+revComp(str(hash2kmer(i,k))))
        bottom.annotate((str(n+1)+". "+str(hash2kmer(i,k))), (xaxis[last], yaxis2c[last]), (xaxis[last]+0.3, yposb[-(n+2)]), size=10, fontname='monospace', weight='bold', color=colourslist[n], arrowprops=dict(color=colourslist[n], shrink=0.05, width=0.05, headwidth=0.4))

    for p, i in enumerate(top6s[::2]):
        dp = (p+3)
        top.text(s=(str(top6labels[p])), x=(numofruns-0.5), y=(ypost[-dp]), size=14, fontname='monospace', color=colourslist[p])

    bar.bar(xaxis, makeyaxis3(k))

    #plt.show()
    plt.savefig('figures/kmerfreq_'+str(identifier)+'_'+str(k), dpi=600)
    plt.close()

def multigrapher(mink, maxk):
    for k in range(mink, maxk+1):
        grapher(k)

multigrapher(mink, maxk)
