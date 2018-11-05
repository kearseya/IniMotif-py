from itertools import groupby
from collections import Counter
import numpy as np
import sys
import os
from Bio.Seq import Seq



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


FileName = "simplefastafile"
runnum = 1
revcompwanted = True

#FileName = input("Fasta File Name:")
#runnum = input("Run number:")

runlists = ["NA",{},{},{},{},{},{},{},{},{},{}]


def CreateKmerList(FileName, runnum, k):
    #runlists.append(runnum)
    #runlists[runnum] = {}
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
