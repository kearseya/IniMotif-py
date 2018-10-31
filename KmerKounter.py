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


#runnum = input("RunNum:")
#filename = input("FileName:")
runnum = 1
filename = "simplefastafile"

frunnumdicts = {1 : "fseqdict1", 2 : "fseqdict2", 3 : "fseqdict3"}
rrunnumdicts = {1 : "rseqdict1", 2 : "rseqdict2", 3 : "rseqdict3"}


fastaFileName = open(filename, "r")
#print(str(fastaFileName.readlines()[-2].strip()[1:]))
#print(kmer2hash(str(fastaFileName.readlines()[-1].strip())))
FileName = fastaFileName

fastaFileName = open(filename, "r")
lastbarcode = fastaFileName.read().splitlines()[-2].strip()[1:]
#print(lastbarcode)



def CreateSeqDict(FileName, runnum):
    fseqdict = frunnumdicts.get(runnum)
    fseqdict = {}
    fastaFileName = open(filename, "r")
    FileName = fastaFileName
    for line in FileName:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            barcode = line[1:]
            if barcode not in fseqdict:
                fseqdict[(barcode)] = []
            continue
        fseq = line
        fseqdict[(barcode)].append(kmer2hash(fseq))
        fastaFileName = open(filename, "r")
        if lastbarcode in list(fseqdict.keys()):
            return fseqdict

fseqdict = CreateSeqDict(FileName, runnum)
print(fseqdict)



def CreateRevSeqDict(FileName, runnum):
    rrunnumdicts.get(runnum)
    rseqdict = {}
    fastaFileName = open(filename, "r")
    FileName = fastaFileName
    for rline in fastaFileName:
        rline = rline.strip()
        if not rline:
            continue
        if rline.startswith(">"):
            rbarcode = rline[1:]
            if rbarcode not in rseqdict:
                rseqdict[(rbarcode+"rev")] = []
            continue
        rseq = rline
        rseqdict[(rbarcode+"rev")].append(kmer2hash(revComp(rseq)))
        fastaFileName = open(filename, "r")
        if lastbarcode + "rev" in list(rseqdict.keys()):
            return rseqdict

rseqdict = CreateRevSeqDict(FileName, runnum)


aseqdict = {**fseqdict,**rseqdict}
print(aseqdict)

aseqvals = aseqdict.values()
aseqvalslist = list(aseqvals)
aseqvalsliststr = [str(i) for i in aseqvalslist]

seqlist = [int(i[1:-1]) for i in aseqvalsliststr]
seqlistcount = Counter(seqlist)

seqcountdict = {**seqlistcount}
print(seqcountdict)



seqs = list(seqcountdict.keys())
print(seqs)

length = 6

def kmercount(l):
    kmerdict = {}
    for i in seqs:
        fastaseq = hash2kmer(i,length)
        for x in range(-1,((length+1)-l)):
            kmers = fastaseq[x:x+l]
            if i not in kmerdict:
                kmerdict[i] = []
                continue
            kmerdict[i].append(kmer2hash(kmers))
    return kmerdict

kmerdict = kmercount(4)
print(kmerdict)
