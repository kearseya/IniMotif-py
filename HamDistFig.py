from test import kmercount
from test import numofruns
import matplotlib.pyplot as plt
import numpy as np
import math

k = int(input("K:"))


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


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def hammer(k):
    khams = []
    for _ in range(numofruns + 1):
        khams.append({})
    for j in range(1,len(kmercount)):
        khams[j] = {k: {} for k in range(k+1)}
        hconsensus = max(kmercount[j][k], key=lambda key: kmercount[j][k][key])
        consensus = hash2kmer(hconsensus, k)
        for x in list(kmercount[j][k].keys()):
            values = hash2kmer(x,k)
            ham = hamming_distance(consensus, values)
            khams[j][ham].update({x:kmercount[j][k][x]})
    return khams

"""
def hammer(k):
    khams = []
    for _ in range(numofruns + 1):
        khams.append({})
    for j in range(1,len(kmercount)):
        khams[j] = {k: {} for k in range(k+1)}
        for i in kmercount[j]:
            hconsensus = max(kmercount[j][i], key=lambda key: kmercount[j][i][key])
            consensus = hash2kmer(hconsensus, i)
            for x in list(kmercount[j][i].keys()):
                values = hash2kmer(x,i)
                ham = hamming_distance(consensus, values)
                khams[j][ham].update({x:kmercount[j][i][x]})
    return khams
"""

khams = hammer(k)



def xaxismaker(run, p):
    x = []
    size = len(khams[run][p])
    for _ in range(size):
        x.append(p)
    return x


def yaxismaker(run, x):
    y = []
    for i in khams[run][x]:
        y.append(khams[run][x][i])
    return y


def top8(k):
    global top8
    top8 = []
    l = (len(kmercount)-1)
    keys = list(kmercount[l][k].keys())
    for i in range(8):
        top8.append(keys[i])
    return top8

#print(khams)

def scatter(run, k):
    top8(k)
    xpos = list(khams[run].keys())
    for x in khams[run]:
        xaxis = xaxismaker(run, x)
        yaxis = yaxismaker(run, x)
        plt.scatter(xaxis, yaxis)

    for x in khams[run]:
        for i in khams[run][x]:
            if i in top8:
                plt.annotate((hash2kmer(i,k) + ' ' + str(khams[run][x][i])), (xpos[x], khams[run][x][i]))


    plt.xlabel("Hamming distance")
    plt.ylabel("Kmer count")
    plt.show()


scatter(1, 6)
