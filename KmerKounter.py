from itertools import groupby
from collections import Counter
import numpy as np
import sys
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


fastaFileName = open("simplefastafile", "r")


fseqdict = {}
for line in fastaFileName:
    line = line.strip()
    #if not line:
            #continue
            #if line.startswith(">"):
                #barcode = line[1:]
            #if barcode not in fseqdict:
                #fseqdict[barcode] = []
            #continue
            #fseq = line
            #fseqdict[barcode].append(fseq)
    if not line:
        continue
    if line.startswith(">"):
        barcode = line[1:]
        if barcode not in fseqdict:
            fseqdict[(barcode)] = []
        continue
    fseq = line
    fseqdict[(barcode)].append(kmer2hash(fseq))


#print(fddict)


rfastaFileName = open("simplefastafile", "r")
rseqdict = {}
for rline in rfastaFileName:
    rline = rline.strip()
    #if not rline:
            #continue
            #if rline.startswith(">"):
                #rbarcode = rline[1:]
            #if rbarcode not in rseqdict:
                #rseqdict[rbarcode] = []
            #continue
            #rseq = rline
            #rseqdict[rbarcode].append(revComp(rseq))
    if not rline:
        continue
    if rline.startswith(">"):
        rbarcode = rline[1:]
        if rbarcode not in rseqdict:
            rseqdict[(rbarcode+"rev")] = []
        continue
    rseq = rline
    rseqdict[(rbarcode+"rev")].append(kmer2hash(revComp(rseq)))

#print(rseqdict)



aseqdict = {**fseqdict,**rseqdict}
print(addict)

a = addict.values()
b = list(a)
c = [str(i) for i in b]
#print(c)
d = Counter(c)
print(d)
