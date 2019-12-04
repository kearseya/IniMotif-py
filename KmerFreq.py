from KmerKounter import identifier

from KmerKounter import kmercount
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import kmer2hash
from KmerKounter import hash2kmer
from KmerKounter import revnuc
from KmerKounter import revComp
from KmerKounter import revcompwanted
from KmerKounter import startround
from KmerKounter import totaldict

from KmerKounter import numofruns
from HamDistFig import top6all

import matplotlib.pyplot as plt
from matplotlib.transforms import TransformedBbox, Bbox
import numpy as np
import math



def makexaxis():
    xaxis = ([])
    for i in range((startround-1), numofruns+(startround-1)):
        xaxis.append(i)
    return xaxis

#print(top6all)



def colours1():
    global colours
    colours = {}
    for k in range(mink, maxk+1):
        top6s = top6all[numofruns][k]
        #consensusnum = kmercount[numofruns][k][top6s[0]]
        colours[k] = list(kmercount[numofruns][k].keys())[6:66]
        #colours[k] = [j[0] for j in colours[k]]
        for i in top6s:
            if i in colours[k]:
                colours[k].remove(i)
    return colours

colours1()



def makeyaxis1a(i, k):
    yaxis1a = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        if i not in top6s or colours[k]:
            try:
                num =  kmercount[l][k][i]
            except:
                num = 0
            if i != rkmer:
                try:
                    rnum = kmercount[l][k][rkmer]
                except:
                    rnum = 0
            else:
                rnum = 0
            f = (num+rnum)/total
            freq = (f/(1-f))
            yaxis1a.append(math.log10(freq))
    return yaxis1a



def makeyaxis1b(i, k):
    yaxis1b = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        if i in colours[k]:
            try:
                num =  kmercount[l][k][i]
            except:
                num = 0
            if i != rkmer:
                try:
                    rnum = kmercount[l][k][rkmer]
                except:
                    rnum = 0
            else:
                rnum = 0
            f = (num+rnum)/total
            freq = (f/(1-f))
            yaxis1b.append(math.log10(freq))
    return yaxis1b


def makeyaxis1c(i, k):
    yaxis1c = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        if i in top6s:
            try:
                num =  kmercount[l][k][i]
            except:
                num = 0
            if i != rkmer:
                try:
                    rnum = kmercount[l][k][rkmer]
                except:
                    rnum = 0
            else:
                rnum = 0
            f = (num+rnum)/total
            freq = (f/(1-f))
            if freq == 0:
                freq = 1
            yaxis1c.append(math.log10(freq))
    return yaxis1c



def makeyaxis2a(i, k):
    yaxis2a = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        if i not in top6s or colours[k]:
            try:
                num =  kmercount[l][k][i]
            except:
                num = 0
            try:
                rnum = kmercount[l][k][rkmer]
            except:
                rnum = 0
            yaxis2a.append((num+rnum)/total)
    return yaxis2a


def makeyaxis2b(i, k):
    yaxis2b = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        if i in colours[k]:
            try:
                num =  kmercount[l][k][i]
            except:
                num = 0
            try:
                rnum = kmercount[l][k][rkmer]
            except:
                rnum = 0
            yaxis2b.append((num+rnum)/total)
    return yaxis2b


def makeyaxis2c(i, k):
    yaxis2c = ([])
    rkmer = kmer2hash(revComp(hash2kmer(i,k)))
    top6s = top6all[numofruns][k]
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        if i in top6s:
            try:
                num =  kmercount[l][k][i]
            except:
                num = 0
            try:
                rnum = kmercount[l][k][rkmer]
            except:
                rnum = 0
            yaxis2c.append((num+rnum)/total)
    return yaxis2c



def makeyaxis3(k):
    yaxis3 = ([])
    for l in range(1, numofruns+1):
        total = totaldict[l][k]
        yaxis3.append(total)
    return yaxis3



def grapher(k):

    xaxis = makexaxis()
    last = (len(xaxis)-1)
    top6s = top6all[numofruns][k]

    fig = plt.figure(figsize=(10,10))
    grid = plt.GridSpec(2,3,wspace=0.4,hspace=0.3)

    top = fig.add_subplot(grid[:-1,:])
    top.set_xlabel("SELEX round")
    top.set_ylabel("log(f/(1-f))")
    top.set_title("Kmer frequency"+', K: '+str(k))
    top.set_xlim([((startround-1)-0.5), (numofruns+startround)])
    top.set_xticks(np.linspace((startround-1)-0.5, (numofruns+startround), num=((2*(numofruns+1))+2), endpoint=True))

    top.spines['right'].set_visible(False)
    top.spines['top'].set_visible(False)

    bottom = fig.add_subplot(grid[-1,:-1])
    bottom.set_xlabel("SELEX round")
    bottom.set_ylabel("f = (kmer/total)")
    bottom.set_title("Kmer frequency")
    bottom.set_xlim([((startround-1)-0.5), (numofruns+(startround-1))])
    bottom.set_xticks(range((startround-1),(numofruns+1)))

    bottom.spines['right'].set_visible(False)
    bottom.spines['top'].set_visible(False)

    bar = fig.add_subplot(grid[-1,-1:])
    bar.set_xlabel("SELEX round")
    bar.set_ylabel("Total kmers")
    bar.set_title("Kmer total distribution")
    bar.set_xticks(makexaxis())

    colourslist = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']

    if revcompwanted == False:
        for l in range(1, (numofruns+1)):
            hashkeys = list(kmercount[l][k].keys())
            for x in range(0, 66):
                i = hashkeys[x]
                num = kmercount[l][k][i]
                rkmer = kmer2hash(revComp(hash2kmer(i,k)))
                try:
                    rnum = kmercount[l][k][rkmer]
                    if num >= rnum:
                        if i in top6s:
                            yaxis1c = makeyaxis1c(i, k)
                            top.plot(xaxis, yaxis1c, color=colourslist[(top6s.index(i)//2)], linewidth=2, marker="s", markevery=None, zorder=(66-x))
                        if i in colours[k]:
                            yaxis1b = makeyaxis1b(i, k)
                            top.plot(xaxis, yaxis1b, linewidth=1,  marker=".", markevery=None, zorder=(66-x))
                except:
                    continue
            if k <= 5:
                step = 1
            if k > 5:
                step = int(round(len(hashkeys)/1000, 0))
            for x in range(66, len(hashkeys), step):
                i = hashkeys[x]
                num = kmercount[l][k][i]
                rkmer = kmer2hash(revComp(hash2kmer(i,k)))
                try:
                    rnum = kmercount[l][k][rkmer]
                    if num >= rnum:
                        yaxis1a = makeyaxis1a(i, k)
                        top.plot(xaxis, yaxis1a, color = '0.75', linestyle='--', linewidth=0.5, marker="x", markevery=None, alpha=0.5, zorder=0)
                except:
                    continue


        for l in range(1, (numofruns+1)):
            hashkeys = list(kmercount[l][k].keys())
            for x in range(0, 66):
                i = hashkeys[x]
                num = kmercount[l][k][i]
                rkmer = kmer2hash(revComp(hash2kmer(i,k)))
                try:
                    rnum = kmercount[l][k][rkmer]
                    if num > rnum:
                        if i in top6s:
                            yaxis2c = makeyaxis2c(i, k)
                            bottom.plot(xaxis, yaxis2c, color=colourslist[(top6s.index(i)//2)], linewidth=2, marker="s", markevery=None, zorder=(66-x))
                        if i in colours[k]:
                            yaxis2b = makeyaxis2b(i, k)
                            bottom.plot(xaxis, yaxis2b, linewidth=1,  marker=".", markevery=None, zorder=66-x)
                except:
                    continue
            if k <= 5:
                step = 1
            if k > 5:
                step = int(round(len(hashkeys)/1000, 0))
            for x in range(66, len(hashkeys), step):
                i = hashkeys[x]
                num = kmercount[l][k][i]
                rkmer = kmer2hash(revComp(hash2kmer(i,k)))
                try:
                    rnum = kmercount[l][k][rkmer]
                    if num > rnum:
                        yaxis2a = makeyaxis2a(i, k)
                        bottom.plot(xaxis, yaxis2a, color = '0.75', linestyle='--', linewidth=0.5, marker="x", markevery=None, alpha=0.8, zorder=0)
                except:
                    continue


    if revcompwanted == True:
        for l in range(1, (numofruns+1)):
            alreadytop = []
            for i in kmercount[l][k]:
                num = kmercount[l][k][i]
                rkmer = kmer2hash(revComp(hash2kmer(i,k)))
                try:
                    rnum = kmercount[l][k][rkmer]
                    if num >= rnum:
                        alreadytop.append(rkmer)
                        if i not in alreadytop:
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


        for l in range(1, numofruns+1):
            alreadybottom = []
            for i in kmercount[l][k]:
                num = kmercount[l][k][i]
                rkmer = kmer2hash(revComp(hash2kmer(i,k)))
                try:
                    rnum = kmercount[l][k][rkmer]
                    if num >= rnum:
                        alreadybottom.append(rkmer)
                        if i not in alreadybottom:
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
        top.annotate((str(n+1)), (xaxis[last], yaxis1c[last]), (xaxis[last]+0.2, ypost[-(n+2)]), size=10, fontname='monospace', weight='bold', arrowprops=dict(color=colourslist[n], shrink=0.05, width=0.05, headwidth=0.4), color=colourslist[n])
        top6labels.append(str(n+1)+". "+str(hash2kmer(i,k))+" / "+revComp(str(hash2kmer(i,k))))
        bottom.annotate((str(n+1)+". "+str(hash2kmer(i,k))), (xaxis[last], yaxis2c[last]), (xaxis[last]+0.3, yposb[-(n+2)]), size=10, fontname='monospace', weight='bold', color=colourslist[n], arrowprops=dict(color=colourslist[n], shrink=0.05, width=0.05, headwidth=0.4))

    for p, i in enumerate(top6s[::2]):
        dp = (p+3)
        top.text(s=(str(top6labels[p])), x=(numofruns-0.5+(startround-1)), y=(ypost[-dp]), size=14, fontname='monospace', color=colourslist[p])

    bar.bar(xaxis, makeyaxis3(k))

    plt.savefig("figures/"+str(identifier)+"/kmer_frequency/kmerfreq_"+str(identifier)+"_"+str(k), dpi=600)
    plt.close()

def multigrapher(mink, maxk):
    for k in range(mink, maxk+1):
        grapher(k)

multigrapher(mink, maxk)
