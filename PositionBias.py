from KmerKounter import identifier

import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
from collections import Counter
from collections import OrderedDict
from Bio import SeqIO
import os

from seaborn import kdeplot

from KmerKounter import numofruns
from KmerKounter import kmercount
from KmerKounter import barcodechecker
from KmerKounter import lvalues
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import filenames
from KmerKounter import ufilenames
from KmerKounter import revcompwanted
from KmerKounter import startround
from KmerKounter import multiround
from KmerKounter import knownbarcode
if multiround == True:
    from KmerKounter import files
    from KmerKounter import datadir
if knownbarcode == True:
    from KmerKounter import barcodes5
    from KmerKounter import barfiveslice

from KmerKounter import barcodeprimers53

from KmerKounter import inputlist

from KmerKounter import top6all

from KmerKounter import extype
from KmerKounter import nobarcode

from WebLogoMod import nmotifs
from WebLogoMod import removelist



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


TSeqNums = {}
LSeqNums = {}
#Barcodevalues = {}
#Primervalues = {}
numoftfbsseq = {}

def seqcountdictinit():
    for x in range(1, numofruns+1):
        TSeqNums[x] = 0
        LSeqNums[x] = 0
        numoftfbsseq[x] = {}
        """
        if knownbarcode == True:
            Barcodevalues[x] = barcodeprimers53[x][0]
        if knownbarcode == False:
            Barcodevalues[x] = barcodechecker(filenames[x])
        """
        for k in range(mink, maxk+1):
            numoftfbsseq[x][k] = 0

seqcountdictinit()

"""
def TSeqCounter(FileName):
    TSeqNum = 0
    fastaFileName = open(FileName, "r")
    line = fastaFileName.readline()
    if line.startswith(">"):
        devision = 2
    if line.startswith("@"):
        devision = 4
    for linenum, lines in enumerate(fastaFileName):
        pass
    TSeqNum = int((linenum+2)//devision)
    return TSeqNum

def LSeqCounter(FileName):
    LSeqNum = 0
    l = lvalues[(ufilenames[FileName])]
    fastaFileName = open(FileName, "r")
    firstline = fastaFileName.readline()
    if firstline.startswith(">"):
        filetpye = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        if len(line) == l and "N" not in line:
            LSeqNum += 1
        else:
            continue
    return LSeqNum

def allTLSeqCounter():
    for i in range(1, numofruns+1):
            x = TSeqCounter(filenames[i])
            y = LSeqCounter(filenames[i])
            b = barcodechecker(filenames[i])
            TSeqNums.update({i:x})
            LSeqNums.update({i:y})
            Barcodevalues.update({i:b})

allTLSeqCounter()
"""

poslist = []
rposlist = []

hamminglist2 = []

consensuslist = []


def listinit():
    for i in range(0,numofruns+1):
        poslist.append({})
        rposlist.append({})
        consensuslist.append({})
    for j in range(1, numofruns+1):
        if extype == "SELEX":
            l = lvalues[j]
        if knownbarcode == False:
            if nobarcode == True:
                avg5 = 0
                avg3 = 0
            if nobarcode == False:
                avg5 = barcodechecker(filenames[j])
                avg3 = avg5
        if knownbarcode == True:
            avg5 = len(barcodeprimers53[j][0])
            avg3 = len(barcodeprimers53[j][1])
        #mink = int(inputlist[((j-1)*5)+7])
        #maxk = int(inputlist[((j-1)*5)+8])
        for k in range(mink, maxk+1):
            poslist[j][k] = {}
            rposlist[j][k] = {}
            consensuslist[j][k] = {}
            for n in range(1, nmotifs+1):
                poslist[j][k][n] = {}
                rposlist[j][k][n] = {}
                consensuslist[j][k][n] = {}
                if extype == "SELEX":
                    for x in range(0,(l+1-k-(avg5+avg3))):
                        poslist[j][k][n].update({x:0})
                        rposlist[j][k][n].update({x:0})
                else:
                    poslist[j][k][n] = []
                    rposlist[j][k][n] = []

listinit()


def hamdictinit(numofruns):
    for _ in range(numofruns + 1):
        hamminglist2.append({})

hamdictinit(numofruns)



def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer():
    for runnum in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            hamminglist2[runnum][k] = {}
            for n in range(1, nmotifs+1):
                hamminglist2[runnum][k][n] = set()
                #for i in kmercount[runnum]:
                    #hconsensus = (list(kmercount[runnum][i].keys())[0])
                if n == 1:
                    hconsensus = max(kmercount[runnum][k], key=lambda key: kmercount[runnum][k][key])
                    consensuslist[runnum][k][n] = hash2kmer(hconsensus, k)
                else:
                    temptop6 = top6all[runnum][k].copy()
                    for j in temptop6:
                        if j in removelist[runnum][k][n]:
                            temptop6.remove(j)
                            if kmer2hash(revComp(hash2kmer(j, k))) in temptop6:
                                temptop6.remove(kmer2hash(revComp(hash2kmer(j, k))))
                    hconsensus = max(temptop6, key=lambda key: kmercount[runnum][k][key])
                    consensuslist[runnum][k][n] = hash2kmer(hconsensus, k)
                consensus = hash2kmer(hconsensus, k)
                for x in list(kmercount[runnum][k].keys()):
                    values = hash2kmer(x,k)
                    rvalues = revComp(values)
                    ham = hamming_distance(consensus, values)
                    rham = hamming_distance(consensus, rvalues)
                    if ham <= 2 and x not in removelist[runnum][k][n]:
                        if kmer2hash(values) not in hamminglist2[runnum][k][n]:
                            hamminglist2[runnum][k][n].add(kmer2hash(values))
                    if rham <= 2 and x not in removelist[runnum][k][n]:
                        if kmer2hash(rvalues) not in hamminglist2[runnum][k][n]:
                            hamminglist2[runnum][k][n].add(kmer2hash(rvalues))

listhammer()

"""
    if revcompwanted == True:
        #for i in kmercount[runnum]:
            #hconsensus = (list(kmercount[runnum][i].keys())[0])
        hconsensus = max(kmercount[runnum][k], key=lambda key: kmercount[runnum][k][key])
        consensus = hash2kmer(hconsensus, k)
        rconsensus = revComp(consensus)
        for x in list(kmercount[runnum][k].keys()):
            values = hash2kmer(x,k)
            rvalues = revComp(values)
            ham = hamming_distance(consensus, values)
            rham = hamming_distance(consensus, rvalues)
            if ham <= 2 or rham <= 2:
                if kmer2hash(rvalues) not in hamminglist2[runnum][k]:
                    hamminglist2[runnum][k].append(kmer2hash(rvalues))
"""


"""
def multilisthammer(numofruns):
    for x in range(1, numofruns+1):
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for i in range(mink, maxk+1):
            listhammer(x, i)

multilisthammer(numofruns)
"""



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



def CreatePosListSELEX(FileName, k, runnum):
    fastaFileName = open(FileName, "r")
    if knownbarcode == False:
        if nobarcode == True:
            avg5 = 0
            avg3 = 0
        if nobarcode == False:
            avg5 = barcodechecker(FileName)
            avg3 = avg5
    if knownbarcode == True:
        avg5 = len(barcodeprimers53[runnum][0])
        avg3 = len(barcodeprimers53[runnum][1])
    combinations = 4**k
    firstline = fastaFileName.readline()
    firstline = firstline.strip()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        c = 0
        if len(line) == lvalues[runnum] and "N" not in line:
            LSeqNums[runnum] += 1
            TSeqNums[runnum] += 1
            done = {}
            for n in range(1, nmotifs+1):
                done[n] = set()
                for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                    kmers = str(line[x+avg5:x+k+avg5])
                    if kmers not in done[n]:
                        hkmers = kmer2hash(kmers)
                        if hkmers < combinations:
                            if hkmers in hamminglist2[runnum][k][n]:
                                poslist[runnum][k][n][x] += 1
                                done[n].add(kmers)
                                if c == 0:
                                    numoftfbsseq[runnum][k] += 1
                                    c += 1
                                #if revcompwanted == True:
                                    #rkmers = revComp(kmers)
                                    #if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                                        #hrkmers = kmer2hash(rkmers)
                                        #if hrkmers < combinations:
                                            #if hrkmers in hamminglist2[runnum][k][n]:
                                                #rposlist[runnum][k][n][x].append(hrkmers)
        else:
            TSeqNums[runnum] += 1
            continue

        if revcompwanted == True:
            if len(line) == lvalues[runnum] and "N" not in line:
                rdone = {}
                for n in range(1, nmotifs+1):
                    rdone[n] = set()
                    for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                        #kmers = str(line[x+avg5:x+k+avg5])
                        rkmers = revComp(line[x+avg5:x+k+avg5])
                        if rkmers not in rdone[n]:
                            hrkmers = kmer2hash(rkmers)
                            if hrkmers in hamminglist2[runnum][k][n]:
                                rposlist[runnum][k][n][x] += 1
                                rdone[n].add(rkmers)



def CreatePosListNORMAL(FileName, k, runnum):
    fastaFileName = open(FileName, "r")
    if knownbarcode == False:
        if nobarcode == True:
            avg5 = 0
            avg3 = 0
        if nobarcode == False:
            avg5 = barcodechecker(FileName)
            avg3 = avg5
    if knownbarcode == True:
        avg5 = len(barcodeprimers53[runnum][0])
        avg3 = len(barcodeprimers53[runnum][1])
    combinations = 4**k
    firstline = fastaFileName.readline()
    firstline = firstline.strip()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        lenline = len(line)
        c = 0
        if "N" not in line:
            LSeqNums[runnum] += 1
            TSeqNums[runnum] += 1
            done = {}
            for n in range(1, nmotifs+1):
                done[n] = set()
                for x in range(0,((lenline+1)-k)-(avg5+avg3)):
                    kmers = str(line[x+avg5:x+k+avg5])
                    if kmers not in done[n]:
                        hkmers = kmer2hash(kmers)
                        if hkmers < combinations:
                            if hkmers in hamminglist2[runnum][k][n]:
                                poslist[runnum][k][n].append(x/lenline)
                                done[n].add(kmers)
                                if c == 0:
                                    numoftfbsseq[runnum][k] += 1
                                    c += 1
                                #if revcompwanted == True:
                                    #rkmers = revComp(kmers)
                                    #if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                                        #hrkmers = kmer2hash(rkmers)
                                        #if hrkmers < combinations:
                                            #if hrkmers in hamminglist2[runnum][k][n]:
                                                #rposlist[runnum][k][n][x].append(hrkmers)
        else:
            TSeqNums[runnum] += 1
            continue

        if revcompwanted == True:
            if "N" not in line:
                rdone = {}
                for n in range(1, nmotifs+1):
                    rdone[n] = set()
                    for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                        #kmers = str(line[x+avg5:x+k+avg5])
                        rkmers = revComp(line[x+avg5:x+k+avg5])
                        if rkmers not in rdone[n]:
                            hrkmers = kmer2hash(rkmers)
                            if hrkmers in hamminglist2[runnum][k][n]:
                                rposlist[runnum][k][n].append(x/lenline)
                                rdone[n].add(rkmers)



def multiPosList(numofruns):
    for x in range(1, numofruns+1):
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for i in range(mink, maxk+1):
            if extype == "SELEX":
                CreatePosListSELEX(filenames[x],i,x)
            else:
                CreatePosListNORMAL(filenames[x],i,x)


def CreatePosListmulti(FileName, k):
    fastaFileName = open(FileName, "r")
    combinations = 4**k
    firstline = fastaFileName.readline()
    firstline = firstline.strip()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        try:
            runnum = barcodes5[line[:barfiveslice]]
            l = lvalues[runnum]
        except:
            continue
        c = 0
        if len(line) == l and "N" not in line:
            LSeqNums[runnum] += 1
            TSeqNums[runnum] += 1
            done = {}
            for n in range(1, nmotifs+1):
                done[n] = set()
                for x in range(0,(l+1-k-len(barcodeprimers53[runnum][0])-len(barcodeprimers53[runnum][1]))):
                    kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                    if kmers not in done[n]:
                        hkmers = kmer2hash(kmers)
                        if hkmers < combinations:
                            if hkmers in hamminglist2[runnum][k][n]:
                                poslist[runnum][k][n][x] += 1
                                done[n].add(kmers)
                                if c == 0:
                                    numoftfbsseq[runnum][k] += 1
                                    c += 1
                                #if revcompwanted == True:
                                    #rkmers = revComp(kmers)
                                    #if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                                        #hrkmers = kmer2hash(rkmers)
                                        #if hrkmers < combinations:
                                            #if hrkmers in hamminglist2[runnum][k][n]:
                                                #rposlist[runnum][k][n][x].append(hrkmers)
        else:
            TSeqNums[runnum] += 1
            continue

        if revcompwanted == True:
            if len(line) == lvalues[runnum] and "N" not in line:
                rdone = {}
                for n in range(1, nmotifs+1):
                    rdone[n] = set()
                    for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                        #kmers = str(line[x+avg5:x+k+avg5])
                        rkmers = revComp(line[x+avg5:x+k+avg5])
                        if rkmers not in rdone[n]:
                            hrkmers = kmer2hash(rkmers)
                            if hrkmers in hamminglist2[runnum][k][n]:
                                rposlist[runnum][k][n][x] += 1
                                rdone[n].add(rkmers)

def multiPosListmulti(numofruns):
    for j in files:
        fullname = os.path.join(datadir, j)
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for i in range(mink, maxk+1):
            CreatePosListmulti(fullname,i)


def PosListMake():
    if multiround == False:
        multiPosList(numofruns)
    if multiround == True:
        multiPosListmulti(numofruns)

PosListMake()


def FindTotal(runnum, k, n):
    total = 0
    for i in poslist[runnum][k][n]:
        total += poslist[runnum][k][n][i]
        if revcompwanted == True:
            total += rposlist[runnum][k][n][i]
    return total



def fseqbias(runnum, k, total, n):
    fseqbias = []
    for i in poslist[runnum][k][n]:
        fseqbias.append(poslist[runnum][k][n][i]/total)
    return fseqbias

def rseqbias(runnum, k, total, n):
    rseqbias = []
    for i in rposlist[runnum][k][n]:
        rseqbias.append(rposlist[runnum][k][n][i]/total)
    return rseqbias



def makexaxis(runnum, k, n):
    xaxis = ([])
    for i in range(1, len(poslist[runnum][k][n])+1):
        xaxis.append(i)
    return xaxis



def rmakexaxis1(width, runnum, k, n):
    fxaxis = ([])
    for i in range(1, len(poslist[runnum][k][n])+1):
        fxaxis.append(i-width/2)
    return fxaxis

def rmakexaxis2(width, runnum, k, n):
    rxaxis = ([])
    for i in range(1, len(poslist[runnum][k][n])+1):
        rxaxis.append(i+width/2)
    return rxaxis



def makeyaxis():
    yaxis = ([])
    for i in range(0, 8):
        yaxis.append(i/10)
    return yaxis



def plotter(runnum, k, n):
    if extype == "SELEX":
        fig, bar = plt.subplots()
        width = 0.4

        xaxis = makexaxis(runnum, k, n)
        yaxis = makeyaxis()

        total = FindTotal(runnum, k, n)

        fseq = fseqbias(runnum, k, total, n)
        rseq = rseqbias(runnum, k, total, n)

        fxaxis = rmakexaxis1(width, runnum, k, n)
        rxaxis = rmakexaxis2(width, runnum, k, n)

        average = (1/len(fseq))
        raverage = (1/(len(fseq)+len(rseq)))

        bar.set_xlabel("TFBS")
        bar.set_ylabel("Frequency")
        bar.set_title("Position Distribution for Run: "+str(runnum+(startround-1))+", K: "+str(k)+", Consensus: "+str(consensuslist[runnum][k][n]))
        bar.set_ylim(0,0.7)
        bar.set_yticks(yaxis)

        if revcompwanted == False:
            bar.set_xticks(xaxis)
            bar.bar(xaxis, fseq)
            bar.axhline(y=average, xmin=0.01, xmax=0.99, linestyle='dashed', color='black')
            bar.text((len(xaxis)-5), 0.65, "Average = "+str(round(average, 4)))

        if revcompwanted == True:
            bar.set_xticks(xaxis)
            bar.bar(fxaxis, fseq, width, label = 'Forward strands')
            bar.bar(rxaxis, rseq, width, label = 'Reverse strands')
            bar.axhline(y=raverage, xmin=0.01, xmax=0.99, linestyle='dashed', color = 'black', linewidth = 0.75)
            bar.text((len(xaxis)-5), 0.55, "Average = "+str(round(raverage, 4)))

            handels=('Forward strands', 'Reverse strands')
            label=('Forward strands', 'Reverse strands')
            bar.legend(loc = 1)

    else:
        fig, kde = plt.subplots()
        kde.set_xlabel("TFBS relative position")
        kde.set_ylabel("Probability density function")
        kde.set_title("TFBS Position Distribution for Run: "+str(runnum+(startround-1))+", K: "+str(k)+", Consensus: "+str(consensuslist[runnum][k][n]))
        kde.set_xlim(0,1)
        if revcompwanted == False:
            plot = kdeplot(poslist[runnum][k][n], shade=True)
        if revcompwanted == True:
            plot1 = kdeplot(poslist[runnum][k][n], shade=True, label='Forward strands')
            plot2 = kdeplot(rposlist[runnum][k][n], shade=True, label='Reverse strands')

    plt.savefig("figures/"+str(identifier)+"/position_bias/pos_"+str(identifier)+"_"+str(runnum+(startround-1))+"_"+str(k)+"_"+str(n), dpi=600)
    plt.close()



def plotrange(numofruns):
    for r in range(1, numofruns+1):
        #mink = int(inputlist[((r-1)*5)+7])
        #maxk = int(inputlist[((r-1)*5)+8])
        for k in range(mink, maxk+1):
            for n in range(1, nmotifs+1):
                plotter(r, k, n)

plotrange(numofruns)



numofkmers = {}
numofuniquekmers = {}
numoftfbs = {}

def numberofukmers():
    for i in range(1, numofruns+1):
        numofkmers.update({i:{}})
        numofuniquekmers.update({i:{}})
        numoftfbs.update({i:{}})
    for i in range(1, numofruns+1):
        #mink = int(inputlist[((i-1)*5)+7])
        #maxk = int(inputlist[((i-1)*5)+8])
        for k in range(mink, maxk+1):
            numofkmers[i].update({k:sum(kmercount[i][k].values())})
            numofuniquekmers[i].update({k:len(kmercount[i][k])})
            numoftfbs[i].update({k:len(hamminglist2[i][k])})

numberofukmers()

"""
numoftfbsseq = {}

def seqwtfbsfinder(FileName, k):
    TFBSSeqNum = 0
    l = lvalues[(ufilenames[FileName])]
    fastaFileName = open(FileName, "r")
    avg = Barcodevalues[(ufilenames[FileName])]
    for line in fastaFileName:
        c = 0
        line = line.strip()
        if line.startswith(">") or line.startswith("@"):
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
                            #if revcompwanted == True:
                                #rkmers = revComp(kmers)
                                #if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                                    #hrkmers = kmer2hash(rkmers)
                                    #if hrkmers < combinations:
                                        #if hrkmers in hamminglist2[(ufilenames[FileName])][k]:
                                            #if c <= 1:
                                                #TFBSSeqNum += 1
                                                #c += 1

        if line.startswith("+"):
            next(fastaFileName)
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
        #mink = int(inputlist[((r-1)*5)+7])
        #maxk = int(inputlist[((r-1)*5)+8])
        for k in range(mink, maxk+1):
            numoftfbsseq[r].update({k:seqwtfbsfinder(filenames[r], k)})

numoftfbsseqs()
"""


def seqbiasfinder():
    seqbias = {}
    for r in range(1, numofruns+1):
        seqbias.update({r:{}})
        #mink = int(inputlist[((r-1)*5)+7])
        #maxk = int(inputlist[((r-1)*5)+8])
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
