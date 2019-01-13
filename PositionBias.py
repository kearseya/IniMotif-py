import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
from collections import Counter
from collections import OrderedDict

from test import numofruns
from test import kmercount

revcompwanted = True

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


def addrun():
    global FileName
    FileName = input("Fasta File Name:")
    global runnum
    runnum = int(input("Run number:"))
    global l
    l = int(input("Read lengths:"))
    global k
    k = int(input("K:"))

addrun()


revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev


def barcodechecker(FileName):
    bar = []
    fastaFileName = open(FileName, "r")
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        for i in range(6):
            start = line[0:i]
            end = revComp(line[len(line)-i : len(line)])
            if start == end:
                x = i
            else:
                break
        bar.append(x)
    favg = (sum(bar)/len(bar))
    avg = round(favg)
    return avg



poslist = []
posdict = []

rposlist = []
rposdict = []
hamminglist2 = []


def listinit(l, k):
    avg = barcodechecker(FileName)
    for _ in range(0,(l+1-k-(2*avg))):
        poslist.append([])
        #posdict.append({})
        rposlist.append([])
        #rposdict.append({})

listinit(l,k)

def hamdictinit(numofruns):
    for _ in range(numofruns + 1):
        hamminglist2.append({})

hamdictinit(numofruns)


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer():
    for i in kmercount[runnum]:
        hamminglist2[runnum][i] = []
        #hconsensus = (list(kmercount[runnum][i].keys())[0])
        hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
        consensus = hash2kmer(hconsensus, i)
        for x in list(kmercount[runnum][i].keys()):
            values = hash2kmer(x,i)
            if hamming_distance(consensus, values) <= 2:
                hamminglist2[runnum][i].append(kmer2hash(values))

    if revcompwanted == True:
        for i in kmercount[runnum]:
            hamminglist2[runnum][i] = []
            #hconsensus = (list(kmercount[runnum][i].keys())[0])
            hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
            consensus = hash2kmer(hconsensus, i)
            rconsensus = revComp(consensus)
            for x in list(kmercount[runnum][i].keys()):
                values = hash2kmer(x,i)
                if hamming_distance(rconsensus, values) <= 2:
                    hamminglist2[runnum][i].append(kmer2hash(values))

listhammer()
print(hamminglist2)




def FindTotal(FileName, k):
    fastaFileName = open(FileName, "r")
    avg = barcodechecker(FileName)
    total = 0
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == l:
            for x in range(0,((len(line)+1)-k)-(2*avg)):
                kmers = str(line[x+avg:x+k+avg])
                if len(kmers) > 0 and line.count(kmers) == 1:
                    hkmers = kmer2hash(kmers)
                    if hkmers in hamminglist2[runnum][k]:
                        total += 1

        if revcompwanted == True:
            if len(line) == l:
                for x in range(0,((len(line)+1)-k)-(2*avg)):
                    rkmers = revComp(line[x+avg:x+k+avg])
                    if len(rkmers) > 0 and (line.count(rkmers) + line.count(kmers)) == 1:
                        hrkmers = kmer2hash(rkmers)
                        if hrkmers in hamminglist2[runnum][k]:
                            total += 1
    return total


total = FindTotal(FileName, k)



def CreatePosList(FileName, k):
    fastaFileName = open(FileName, "r")
    avg = barcodechecker(FileName)
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == l:
            for x in range(0,((len(line)+1)-k)-(2*avg)):
                kmers = str(line[x+avg:x+k+avg])
                if len(kmers) > 0 and line.count(kmers) == 1:
                    hkmers = kmer2hash(kmers)
                    if hkmers in hamminglist2[runnum][k]:
                        poslist[x].append(hkmers)

        if revcompwanted == True:
            if len(line) == l:
                for x in range(0,((len(line)+1)-k)-(2*avg)):
                    rkmers = revComp(line[x+avg:x+k+avg])
                    if len(rkmers) > 0 and (line.count(rkmers) + line.count(kmers)) == 1:
                        hrkmers = kmer2hash(rkmers)
                        if hrkmers in hamminglist2[runnum][k]:
                            rposlist[x].append(hrkmers)

CreatePosList(FileName, k)

"""
def KPosCounter():
    for x in range(0,len(poslist)):
        poslistcount = Counter(poslist[x])
        poslistcount2 = {**poslistcount}
        posdict[x].update(poslistcount2)

def rKPosCounter():
    for x in range(0,len(rposlist)):
        rposlistcount = Counter(rposlist[x])
        rposlistcount2 = {**rposlistcount}
        rposdict[x].update(rposlistcount2)
"""
"""
KPosCounter()
rKPosCounter()
"""

print(rposlist[0])
print(rposlist[1])




def fseqbias():
    fseqbias = []
    for i in poslist:
        fseqbias.append(len(i)/total)
    return fseqbias

def rseqbias():
    rseqbias = []
    for i in rposlist:
        rseqbias.append(len(i)/total)
    return rseqbias


fseq = fseqbias()
rseq = rseqbias()

print(fseq)
print(sum(fseq))
print(rseq)
print(sum(rseq))

print(1/(len(fseq)+len(rseq)))




def makexaxis():
    xaxis = ([])
    for i in range(1, len(poslist)+1):
        xaxis.append(i)
    return xaxis

#xaxis = makexaxis()


def rmakexaxis1(width):
    fxaxis = ([])
    for i in range(1, len(poslist)+1):
        fxaxis.append(i-width/2)
    return fxaxis

def rmakexaxis2(width):
    rxaxis = ([])
    for i in range(1, len(poslist)+1):
        rxaxis.append(i+width/2)
    return rxaxis



def makeyaxis():
    yaxis = ([])
    for i in range(0, 8):
        yaxis.append(i/10)
    return yaxis



def plotter():
    fig, bar = plt.subplots()
    width = 0.4

    xaxis = makexaxis()
    yaxis = makeyaxis()

    fxaxis = rmakexaxis1(width)
    rxaxis = rmakexaxis2(width)

    average = (1/len(fseq))
    raverage = (1/(len(fseq)+len(rseq)))

    bar.set_xlabel("TFBS")
    bar.set_ylabel("Frequency")
    bar.set_title("Position Distribution")
    bar.set_ylim(0,0.7)
    bar.set_yticks(yaxis)

    if revcompwanted == False:
        bar.set_xticks(xaxis)
        bar.bar(xaxis, fseq)
        bar.axhline(y=average, xmin=0.05, xmax=0.95, linestyle='dashed', color='black')

    if revcompwanted == True:
        bar.set_xticks(xaxis)
        bar.bar(fxaxis, fseq, width, label = 'Forward strands')
        bar.bar(rxaxis, rseq, width, label = 'Reverse strands')
        bar.axhline(y=raverage, xmin=0.01, xmax=0.99, linestyle='dashed', color = 'black')



    handels=('Forward strands', 'Reverse strands')
    label=('Forward strands', 'Reverse strands')
    bar.legend(loc = 1)
    plt.show()


plotter()













"""
            try:
                ftotal = kmercount[x+1][k][i]
                #print(ftotal)
            except:
                continue
            #print(ftotal)
            try:
                rtotal = kmercount[x+1][k][rseq]
                #print(rtotal)
            except:
                continue
"""


def TSeqCounter(FileName, k):
    global TSeqNum
    TSeqNum = 0
    fastaFileName = open(FileName, "r")
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        TSeqNum += 1

def LSeqCounter(FileName, k):
    global LSeqNum
    LSeqNum = 0
    fastaFileName = open(FileName, "r")
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == l:
            LSeqNum += 1
    #return SeqNum
"""
TSeqCounter(FileName, k)
LSeqCounter(FileName, k)
print(TSeqNum)
print(LSeqNum)
"""
