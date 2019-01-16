import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
from collections import Counter
from collections import OrderedDict

from KmerKount import numofruns
from KmerKount import kmercount
from KmerKount import barcodechecker
from KmerKount import lvalues
from KmerKount import mink
from KmerKount import maxk
from KmerKount import filenames
from KmerKount import revcompwanted



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

#addrun()


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
    for i in kmercount[runnum]:
        hamminglist2[runnum][k] = []
        #hconsensus = (list(kmercount[runnum][i].keys())[0])
        hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
        consensus = hash2kmer(hconsensus, i)
        for x in list(kmercount[runnum][i].keys()):
            values = hash2kmer(x,i)
            rvalues = revComp(values)
            ham = hamming_distance(consensus, values)
            rham = hamming_distance(consensus, rvalues)
            if ham or rham <= 2:
                hamminglist2[runnum][k].append(kmer2hash(values))

    if revcompwanted == True:
        for i in kmercount[runnum]:
            #hamminglist2[runnum][k] = []
            #hconsensus = (list(kmercount[runnum][i].keys())[0])
            hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
            consensus = hash2kmer(hconsensus, i)
            rconsensus = revComp(consensus)
            for x in list(kmercount[runnum][i].keys()):
                values = hash2kmer(x,i)
                if hamming_distance(rconsensus, values) <= 2:
                    hamminglist2[runnum][k].append(kmer2hash(values))




def multilisthammer(numofruns, mink, maxk):
    for x in range(1, numofruns+1):
        for i in range(mink, maxk+1):
            listhammer(x, i)

multilisthammer(numofruns, mink, maxk)




"""
def FindTotal(FileName, k, runnum):
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
            #total = total*2
            if len(line) == l:
                for x in range(0,((len(line)+1)-k)-(2*avg)):
                    rkmers = revComp(line[x+avg:x+k+avg])
                    if len(rkmers) > 0 and (line.count(rkmers) + line.count(kmers)) == 1:
                        hrkmers = kmer2hash(rkmers)
                        if hrkmers in hamminglist2[runnum][k]:
                            total += 1

    return total
"""



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
        #print(poslist[runnum][k][i])
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
        bar.axhline(y=raverage, xmin=0.01, xmax=0.99, linestyle='dashed', color = 'black')
        bar.text((len(xaxis)-5), 0.55, "Average = "+str(round(raverage, 4)))

        handels=('Forward strands', 'Reverse strands')
        label=('Forward strands', 'Reverse strands')
        bar.legend(loc = 1)

    plt.show()


#plotter()



def plotrange(numofruns, mink, maxk):
    for r in range(1,numofruns+1):
        for k in range(mink, maxk+1):
            plotter(r, k)

plotrange(numofruns, mink, maxk)



































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




"""
for x in range(0, len(poslist)):
    for i in poslist[x]:
        fseq = hash2kmer(i,k)
        rseq = kmer2hash(revComp(fseq))
        fseqcount = posdict[x][i]
    #    print(fseqcount)
        try:
            rseqcount = posdict[x][rseq]
            #print(rseqcount)
        except:
            rseqcount = 0
    #    print(rseqcount)

        #print(rtotal)
        bias = (((fseqcount)-(rseqcount))/((fseqcount)+(rseqcount)))
        #print(bias)
        #print(fseqbias)
        fseqbias[x].append(bias)
        #print(fseqbias)
print(fseqbias)
for x in range(0,len(poslist)):
    fseqbiasval = (sum(fseqbias[x])/len(fseqbias[x]))
    fseqvals.append(fseqbiasval)
return fseqvals
"""
