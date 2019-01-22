from KmerKounter import identifier

import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
from collections import Counter
from collections import OrderedDict

from KmerKounter import numofruns
from KmerKounter import runlists
from KmerKounter import kmercount
from KmerKounter import barcodechecker
from KmerKounter import lvalues
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import filenames
from KmerKounter import ufilenames
from KmerKounter import revcompwanted



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

"""
def addrun():
    global FileName
    FileName = input("Fasta File Name:")
    global runnum
    runnum = int(input("Run number:"))
    global l
    l = int(input("Read lengths:"))
    global k
    k = int(input("K:"))
"""


revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev


poslist = []
rposlist = []

hamminglist2 = []


def listinit(mink, maxk):
    l = lvalues[1]
    for i in range(0,numofruns+1):
        poslist.append({})
        rposlist.append({})
        try:
            avg = barcodechecker(filenames[i])
        except:
            continue
        for k in range(mink, maxk+1):
            poslist[i].update({k:{}})
            rposlist[i].update({k:{}})
            for x in range(0,(l+1-k-(2*avg))):
                poslist[i][k].update({x:[]})
                rposlist[i][k].update({x:[]})


listinit(mink,maxk)


def hamdictinit(numofruns):
    for _ in range(numofruns + 1):
        hamminglist2.append({})

hamdictinit(numofruns)



def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer(runnum, k):
    hamminglist2[runnum][k] = []
    hconsensus = max(kmercount[runnum][k], key=lambda key: kmercount[runnum][k][key])
    consensus = hash2kmer(hconsensus, k)
    for x in list(kmercount[runnum][k].keys()):
        values = hash2kmer(x,k)
        rvalues = revComp(values)
        ham = hamming_distance(consensus, values)
        rham = hamming_distance(consensus, rvalues)
        if ham <= 2:
            if kmer2hash(values) not in hamminglist2[runnum][k]:
                hamminglist2[runnum][k].append(kmer2hash(values))
        if rham <= 2:
            if kmer2hash(rvalues) not in hamminglist2[runnum][k]:
                hamminglist2[runnum][k].append(kmer2hash(rvalues))




def multilisthammer(numofruns, mink, maxk):
    for x in range(1, numofruns+1):
        for i in range(mink, maxk+1):
            listhammer(x, i)

multilisthammer(numofruns, mink, maxk)



def CreatePosList(FileName, k, runnum):
    fastaFileName = open(FileName, "r")
    avg = barcodechecker(FileName)
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == lvalues[runnum] and "N" not in line:
            for x in range(0,((len(line)+1)-k)-(2*avg)):
                kmers = str(line[x+avg:x+k+avg])
                if len(kmers) > 0 and line.count(kmers) == 1:
                    hkmers = kmer2hash(kmers)
                    if hkmers in hamminglist2[runnum][k]:
                        poslist[runnum][k][x].append(hkmers)

        if revcompwanted == True:
            if len(line) == lvalues[runnum] and "N" not in line:
                for x in range(0,((len(line)+1)-k)-(2*avg)):
                    kmers = str(line[x+avg:x+k+avg])
                    rkmers = revComp(line[x+avg:x+k+avg])
                    if len(rkmers) > 0 and (line.count(rkmers) + line.count(kmers)) == 1:
                        hrkmers = kmer2hash(rkmers)
                        if hrkmers in hamminglist2[runnum][k]:
                            rposlist[runnum][k][x].append(hrkmers)


def multiPosList(numofruns, mink, maxk):
    for x in range(1, numofruns+1):
        for i in range(mink, maxk+1):
            CreatePosList(filenames[x],i,x)

multiPosList(numofruns, mink, maxk)




def FindTotal(runnum, k):
    total = 0
    for i in poslist[runnum][k]:
        total += len(poslist[runnum][k][i])
        if revcompwanted == True:
            total += len(rposlist[runnum][k][i])
    return total



def fseqbias(runnum, k, total):
    fseqbias = []
    for i in poslist[runnum][k]:
        fseqbias.append(len(poslist[runnum][k][i])/total)
    return fseqbias

def rseqbias(runnum, k, total):
    rseqbias = []
    for i in rposlist[runnum][k]:
        rseqbias.append(len(rposlist[runnum][k][i])/total)
    return rseqbias



def makexaxis(runnum, k):
    xaxis = ([])
    for i in range(1, len(poslist[runnum][k])+1):
        xaxis.append(i)
    return xaxis



def rmakexaxis1(width, runnum, k):
    fxaxis = ([])
    for i in range(1, len(poslist[runnum][k])+1):
        fxaxis.append(i-width/2)
    return fxaxis

def rmakexaxis2(width, runnum, k):
    rxaxis = ([])
    for i in range(1, len(poslist[runnum][k])+1):
        rxaxis.append(i+width/2)
    return rxaxis



def makeyaxis():
    yaxis = ([])
    for i in range(0, 8):
        yaxis.append(i/10)
    return yaxis



def plotter(runnum, k):
    fig, bar = plt.subplots()
    width = 0.4

    xaxis = makexaxis(runnum, k)
    yaxis = makeyaxis()

    total = FindTotal(runnum, k)

    fseq = fseqbias(runnum, k, total)
    rseq = rseqbias(runnum, k, total)

    fxaxis = rmakexaxis1(width, runnum, k)
    rxaxis = rmakexaxis2(width, runnum, k)

    average = (1/len(fseq))
    raverage = (1/(len(fseq)+len(rseq)))

    bar.set_xlabel("TFBS")
    bar.set_ylabel("Frequency")
    bar.set_title("Position Distribution for Run: "+str(runnum)+", K: "+str(k))
    bar.set_ylim(0,0.7)
    bar.set_yticks(yaxis)

    if revcompwanted == False:
        bar.set_xticks(xaxis)
        bar.bar(xaxis, fseq)
        bar.axhline(y=average, xmin=0.01, xmax=0.99, linestyle='dashed', color='black')
        bar.text((len(xaxis)-0.8), 0.55, "Average = "+str(round(average, 4)))

    if revcompwanted == True:
        bar.set_xticks(xaxis)
        bar.bar(fxaxis, fseq, width, label = 'Forward strands')
        bar.bar(rxaxis, rseq, width, label = 'Reverse strands')
        bar.axhline(y=raverage, xmin=0.01, xmax=0.99, linestyle='dashed', color = 'black', linewidth = 0.75)
        bar.text((len(xaxis)-5), 0.55, "Average = "+str(round(raverage, 4)))

        handels=('Forward strands', 'Reverse strands')
        label=('Forward strands', 'Reverse strands')
        bar.legend(loc = 1)

    plt.savefig('figures/pos_'+str(identifier)+"_"+str(runnum)+"_"+str(k), dpi=600)
    plt.close()



def plotrange(numofruns, mink, maxk):
    for r in range(1,numofruns+1):
        for k in range(mink, maxk+1):
            plotter(r, k)

plotrange(numofruns, mink, maxk)




def TSeqCounter(FileName):
    TSeqNum = 0
    fastaFileName = open(FileName, "r")
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        TSeqNum += 1
    return TSeqNum

def LSeqCounter(FileName):
    LSeqNum = 0
    l = lvalues[(ufilenames[FileName])]
    fastaFileName = open(FileName, "r")
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == l and "N" not in line:
            LSeqNum += 1
    return LSeqNum


TSeqNums = {}
LSeqNums = {}
Barcodevalues = {}


def allTLSeqCounter():
    for i in range(1, numofruns+1):
        x = TSeqCounter(filenames[i])
        y = LSeqCounter(filenames[i])
        b = barcodechecker(filenames[i])
        TSeqNums.update({i:x})
        LSeqNums.update({i:y})
        Barcodevalues.update({i:b})

allTLSeqCounter()

numofkmers = {}
numofuniquekmers = {}
numoftfbs = {}

def numberofukmers():
    for i in range(1, numofruns+1):
        numofkmers.update({i:{}})
        numofuniquekmers.update({i:{}})
        numoftfbs.update({i:{}})
    for i in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            numofkmers[i].update({k:len(runlists[i][k])})
            numofuniquekmers[i].update({k:len(kmercount[i][k])})
            numoftfbs[i].update({k:len(hamminglist2[i][k])})

numberofukmers()

numoftfbsseq = {}

def seqwtfbsfinder(FileName, k):
    TFBSSeqNum = 0
    l = lvalues[(ufilenames[FileName])]
    fastaFileName = open(FileName, "r")
    avg = barcodechecker(FileName)
    for line in fastaFileName:
        c = 0
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == l and "N" not in line:
            for x in range(0,((len(line)+1)-k)-(2*avg)):
                kmers = str(line[x+avg:x+k+avg])
                if len(kmers) > 0 and line.count(kmers) == 1:
                    hkmer = kmer2hash(kmers)
                    if hkmer in hamminglist2[(ufilenames[FileName])][k]:
                        if c == 0:
                            TFBSSeqNum += 1
                            c += 1

        else:
            continue

        if revcompwanted == True:
            if len(line) == l and "N" not in line:
                for x in range(0,((len(line)+1)-k)-(2*avg)):
                    rkmers = revComp(line[x+avg:x+k+avg])
                    if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) == 1:
                         hkmer = kmer2hash(rkmers)
                         if hkmer in hamminglist2[(ufilenames[FileName])][k]:
                             if c == 0:
                                 TFBSSeqNum += 1
                                 c += 1
            else:
                continue
    return TFBSSeqNum

def numoftfbsseqs():
    for i in range(1, numofruns+1):
        numoftfbsseq.update({i:{}})
    for r in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            numoftfbsseq[r].update({k:seqwtfbsfinder(filenames[r], k)})

numoftfbsseqs()



def seqbiasfinder():
    seqbias = {}
    for r in range(1, numofruns+1):
        seqbias.update({r:{}})
        for k in range(mink, maxk+1):
            seqbias[r].update({k:[]})
            for x in kmercount[r][k]:
                try:
                    forward = kmercount[r][k][x]
                    reverse = kmercount[r][k][kmer2hash(revComp(hash2kmer(x,k)))]
                    seqbiasval = ((forward-reverse)/(forward+reverse))
                    if seqbiasval > 0:
                        seqbias[r][k].append(seqbiasval)
                    if seqbiasval < 0:
                        seqbias[r][k].append(val*(-1))
                except:
                    continue
    return seqbias

seqbias = seqbiasfinder()
