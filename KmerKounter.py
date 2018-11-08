from itertools import groupby
from collections import Counter
from collections import OrderedDict
import numpy as np
#import sys
#import os
#from Bio.Seq import Seq



def kmer2hash(kmer):
    """
    kmer: input sequence
    return: a hash code, 32 bit
    """
    k = len(kmer)
    assert k<16, "kmer should be shorted than 16 bases"
    base = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    kh = np.zeros(1,dtype='uint32')
    kh = base[kmer[0]]
    for tb in kmer[1:]:
        kh = kh<<2
        kh += base[tb]

    return kh

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



revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev


FileName = "NF1-2"
runnum = 1
revcompwanted = False

#FileName = input("Fasta File Name:")
#runnum = input("Run number:")

runlists = ["NA",{},{},{},{},{},{},{},{},{},{}]
kmercount = ["NA",{},{},{},{},{},{},{},{},{},{}]
hamlist = ["NA",{},{},{},{},{},{},{},{},{},{}]
hamdict = ["NA",{},{},{},{},{},{},{},{},{},{}]



def CreateKmerList(FileName, runnum, k):
    runlists[runnum][k] = []
    fastaFileName = open(FileName, "r")
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        for x in range(-1,((len(line)+1)-k)):
            kmers = line[x:x+k]
            if len(kmers) > 0:
                runlists[runnum][k].append(kmer2hash(kmers))

        if revcompwanted == True:
            for x in range(-1,((len(line)+1)-k)):
                rkmers = revComp(line[x:x+k])
                if len(rkmers) > 0:
                    runlists[runnum][k].append(kmer2hash(rkmers))



def RangeKmerList(mink,maxk):
    for i in range(mink,maxk+1):
        CreateKmerList(FileName, runnum, i)


RangeKmerList(4,7)



def KmerCounter():
    for x in list(runlists[runnum]):
        kmercount[runnum][x] = {}
    for i in runlists[runnum]:
        kmerlistcount = Counter(runlists[runnum][i])
        kmerlistcount2 = {**kmerlistcount}
        sort = sorted(kmerlistcount2.items(), key=lambda x: x[1], reverse=True)
        sortdict = OrderedDict(sort)
        sortdict2 = {**sortdict}
        kmercount[runnum][i].update(sortdict2)


KmerCounter()

#print(list(kmercount[runnum][5].keys())[0])



def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer():
    for i in kmercount[runnum]:
        hamlist[runnum][i] = []
        hconsensus = (list(kmercount[runnum][i].keys())[0])
        consensus = hash2kmer(hconsensus, i)
        for x in list(kmercount[runnum][i].keys()):
            values = hash2kmer(x,i)
            if hamming_distance(consensus, values) <= 1:
                hamlist[runnum][i].append(kmer2hash(values))


listhammer()



def dicthammer():
    for i in hamlist[runnum]:
        hamdict[runnum][i] = {}
        test = { z : kmercount[runnum][i][z] for z in hamlist[runnum][i] }
        M = sum(test.values())
        test = { z : (kmercount[runnum][i][z]/M) for z in test }
        hamdict[runnum][i].update(test)


dicthammer()


print(hamdict[runnum])
