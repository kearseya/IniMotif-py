from test import kmercount
#from test import hamlist
import matplotlib.pyplot as plt
import numpy as np
import math


def hash2kmer(hashkey, k):
    """
    hashkey: hash key of kmer, numpy, 'uint32'
    k: length of kmer
    """
    base = np.array('ACGT', 'c')
    arr = np.chararray(k)
    mask = 0b11

    arr[-1]=base[ mask & hashkey]
    for i in range(2,k+1):
        hashkey = (hashkey>>2)
        arr[-i]=base[mask & hashkey]

    return arr.tostring().decode("utf-8")


def makexaxis():
    xaxis = ([])
    for i in range(len(kmercount)-1):
        xaxis.append(i)
    return xaxis

"""
def makexaxis2():
    xaxis = ([])
    for i in range(-1, len(kmercount)+1):
        xaxis.append(i)
    return xaxis
"""

def makeyaxis1(i, k):
    yaxis1 = ([])
    for l in range(1, (len(kmercount))):
        total = sum(kmercount[l][k].values())
        if i in kmercount[l][k]:
            num =  kmercount[l][k][i]
            f = num/total
            freq = (f/(1-f))
            yaxis1.append(math.log10(freq))
    return yaxis1

def makeyaxis2(i, k):
    yaxis2 = ([])
    for l in range(1, (len(kmercount))):
        total = sum(kmercount[l][k].values())
        if i in kmercount[l][k]:
            num =  kmercount[l][k][i]
            f = num/total
            yaxis2.append(f)
    return yaxis2

def makeyaxis3(k):
    yaxis3 = ([])
    for l in range(1, (len(kmercount))):
        total = sum(kmercount[l][k].values())
        yaxis3.append(total)
    return yaxis3


def colours(k):
    global colours
    colours = []
    for l in range(1, len(kmercount)):
        keys = list(kmercount[l][k].keys())
        for i in range(25):
            colours.append(keys[i])
    return colours

def top6(k):
    global top6
    top6 = []
    l = (len(kmercount)-1)
    keys = list(kmercount[l][k].keys())
    for i in range(6):
        top6.append(keys[i])
    return top6


#print(colours(6))
#print(top6(6))


def graphyboy(k):

    xaxis = makexaxis()
    last = (len(xaxis)-1)
    top6(k)
    colours(k)

    #graph = plt.subplots(1,2)
    fig = plt.figure(figsize=(6,6))
    grid = plt.GridSpec(2,3,wspace=0.4,hspace=0.3)

    top = fig.add_subplot(grid[:-1,:])
    top.set_xlabel("SELEX round")
    top.set_ylabel("Kmer count (log(f/(1-f)))")
    top.set_title("Kmer frequency")
    top.set_xticks(makexaxis())

    bottom = fig.add_subplot(grid[-1,:-1])
    bottom.set_xlabel("SELEX round")
    bottom.set_ylabel("Kmer count (kmer/total)")
    bottom.set_title("Kmer frequency")
    bottom.set_xticks(makexaxis())

    bar = fig.add_subplot(grid[-1,-1:])
    bar.set_xlabel("SELEX round")
    bar.set_ylabel("Total kmers")
    bar.set_title("Kmer total distribution")
    bar.set_xticks(makexaxis())


    for l in range(1, len(kmercount)):
        for i in kmercount[l][k]:
            yaxis1 = makeyaxis1(i,k)
            if i in colours:
                try:
                    top.plot(xaxis, yaxis1)
                except:
                    continue
            else:
                try:
                    top.plot(xaxis, yaxis1, color = '0.75')
                except:
                    continue
            if i in top6:
                top.annotate((str((top6.index(i)+1)) + '. ' + hash2kmer(i,k)), (xaxis[last], yaxis1[last]))


    for l in range(1, len(kmercount)):
        for i in kmercount[l][k]:
            yaxis2 = makeyaxis2(i,k)
            if i in colours:
                try:
                    bottom.plot(xaxis, yaxis2)
                except:
                    continue
            else:
                try:
                    bottom.plot(xaxis, yaxis2, color = '0.75')
                except:
                    continue
            if i in top6:
                bottom.annotate((str((top6.index(i)+1)) + '. ' + hash2kmer(i,k)), (xaxis[last], yaxis2[last]))

    bar.bar(xaxis, makeyaxis3(k))

    handles = []
    for i in top6:
        handles.append(str((top6.index(i)+1)) + '. ' + hash2kmer(i,k))
    handles.reverse()
    

    top.legend(handles, bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)


    plt.show()


graphyboy(6)
