from KmerKounter import identifier

import numpy as np
import matplotlib.pyplot as plt
from itertools import groupby
from collections import Counter
from collections import OrderedDict
from Bio import SeqIO
import os

from sklearn.neighbors import KernelDensity
import statistics

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
from WebLogoMod import allowham
from WebLogoMod import removelist
from WebLogoMod import logoform



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


TSeqNums = {}
LSeqNums = {}
numoftfbsseq = {}

def seqcountdictinit():
    for x in range(1, numofruns+1):
        TSeqNums[x] = 0
        LSeqNums[x] = 0
        numoftfbsseq[x] = {}
        for k in range(mink, maxk+1):
            numoftfbsseq[x][k] = 0

seqcountdictinit()


poslist = []
rposlist = []

hamminglist2 = []
rhamminglist2 = []

consensuslist = []
cooccurencelist = []
cooccurenceindexlist = []
gaplist = []
gapstdev = []


def listinit():
    for i in range(0,numofruns+1):
        poslist.append({})
        rposlist.append({})
        consensuslist.append({})
        cooccurencelist.append({})
        cooccurenceindexlist.append({})
        gaplist.append({})
        gapstdev.append({})
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
        for k in range(mink, maxk+1):
            poslist[j][k] = {}
            rposlist[j][k] = {}
            consensuslist[j][k] = {}
            cooccurencelist[j][k] = {}
            cooccurenceindexlist[j][k] = {}
            gaplist[j][k] = {}
            gapstdev[j][k] = {}
            for n in range(1, nmotifs+1):
                poslist[j][k][n] = {}
                rposlist[j][k][n] = {}
                consensuslist[j][k][n] = {}
                cooccurencelist[j][k][n] = {}
                cooccurenceindexlist[j][k][n] = {}
                gaplist[j][k][n] = []
                gapstdev[j][k][n] = {}
                gapstdev[j][k][n]["mean"] = 0
                gapstdev[j][k][n]["stdev"] = 0
                if extype == "SELEX":
                    for x in range(0,(l+1-k-(avg5+avg3))):
                        poslist[j][k][n].update({x:0})
                        rposlist[j][k][n].update({x:0})
                else:
                    poslist[j][k][n] = []
                    rposlist[j][k][n] = []

listinit()

def cooccinit():
    for r in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            for n in range(1, nmotifs+1):
                cooccurencelist[r][k][n]["fr"] = 0
                cooccurencelist[r][k][n]["f"] = 0
                cooccurencelist[r][k][n]["r"] = 0

cooccinit()

def hamdictinit(numofruns):
    for _ in range(numofruns + 1):
        hamminglist2.append({})
        rhamminglist2.append({})

hamdictinit(numofruns)



def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer():
    for runnum in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            hamminglist2[runnum][k] = {}
            rhamminglist2[runnum][k] = {}
            for n in range(1, nmotifs+1):
                hamminglist2[runnum][k][n] = set()
                rhamminglist2[runnum][k][n] = set()
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
                        hamminglist2[runnum][k][n].add(x)
                        if rham <= 2 and x not in removelist[runnum][k][n]:
                            hamminglist2[runnum][k][n].add(kmer2hash(rvalues))
                        if rham > 2 and x not in removelist[runnum][k][n]:
                            rhamminglist2[runnum][k][n].add(kmer2hash(rvalues))

listhammer()



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

    firstline = fastaFileName.readline()
    firstline = firstline.strip()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"

    l = lvalues[runnum]
    for n in range(1, nmotifs+1):
        seq1 = kmer2hash(consensuslist[runnum][k][n])
        seq2 = kmer2hash(revComp(consensuslist[runnum][k][n]))
        for sequence in SeqIO.parse(FileName, filetype):
            line = str(sequence.seq)
            c = 0
            if len(line) == l and "N" not in line:
                LSeqNums[runnum] += 1
                TSeqNums[runnum] += 1
                done = set()
                noisefilter = {"F": [],"R": []}
                for x in range(0,((l+1)-k)-(avg5+avg3)):
                    kmers = kmer2hash(str(line[x+avg5:x+k+avg5]))
                    #if kmers not in done:
                    if kmers in hamminglist2[runnum][k][n]:
                        poslist[runnum][k][n][x] += 1
                        done.add(kmers)
                        if c == 0:
                            numoftfbsseq[runnum][k] += 1
                            c += 1
                    if kmers in rhamminglist2[runnum][k][n]:
                        rposlist[runnum][k][n][x] += 1
                        done.add(kmers)
                        if c == 0:
                            numoftfbsseq[runnum][k] += 1
                            c += 1
                    if kmers == seq1:
                        noisefilter["F"].append(x)
                    if kmers == seq2:
                        noisefilter["R"].append(x)
                    if x == int(l-k)-(avg5+avg3):
                        if seq1 & seq2 in done:
                            numf = len(noisefilter["F"])
                            numr = len(noisefilter["R"])
                            a = 0
                            b = 0
                            mindiff = 1000
                            while (a < numf and b < numr):
                                if (k < abs(noisefilter["F"][a] - noisefilter["R"][b]) < mindiff):
                                    mindiff = abs(noisefilter["F"][a] - noisefilter["R"][b])
                                if (noisefilter["F"][a] < noisefilter["R"][b]):
                                    a += 1
                                else:
                                    b += 1
                            if k < mindiff < l:
                                cooccurencelist[runnum][k][n]["fr"] += 1
                                gaplist[runnum][k][n].append(mindiff-k)
                        elif seq1 in done and seq2 not in done:
                            cooccurencelist[runnum][k][n]["f"] += 1
                        elif seq2 in done and seq1 not in done:
                            cooccurencelist[runnum][k][n]["r"] += 1

        else:
            TSeqNums[runnum] += 1
            continue



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

    firstline = fastaFileName.readline()
    firstline = firstline.strip()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"

    for n in range(1, nmotifs+1):
        seq1 = kmer2hash(consensuslist[runnum][k][n])
        seq2 = kmer2hash(revComp(consensuslist[runnum][k][n]))
        for sequence in SeqIO.parse(FileName, filetype):
            line = str(sequence.seq)
            lenline = len(line)
            c = 0
            LSeqNums[runnum] += 1
            TSeqNums[runnum] += 1
            done = set()
            noisefilter = {"F": [],"R": []}
            for x in range(0,(lenline+1-k)-(avg5+avg3)):
                try:
                    kmers = kmer2hash(str(line[x+avg5:x+k+avg5]))
                    #if kmers not in done:
                    if kmers in hamminglist2[runnum][k][n]:
                        poslist[runnum][k][n].append(x/(lenline-k+1-(avg5+avg3)))
                        done.add(kmers)
                        if c == 0:
                            numoftfbsseq[runnum][k] += 1
                            c += 1
                    if kmers in rhamminglist2[runnum][k][n]:
                        rposlist[runnum][k][n].append((x)/(lenline-k+1-(avg5+avg3)))
                        done.add(kmers)
                        if c == 0:
                            numoftfbsseq[runnum][k] += 1
                            c += 1
                    if kmers == seq1:
                        noisefilter["F"].append(x)
                    if kmers == seq2:
                        noisefilter["R"].append(x)
                except:
                    continue
                if x == (len(line)-k)-(avg5+avg3):
                    if seq1 & seq2 in done:
                        numf = len(noisefilter["F"])
                        numr = len(noisefilter["R"])
                        a = 0
                        b = 0
                        mindiff = 1000
                        while (a < numf and b < numr):
                            if (k < abs(noisefilter["F"][a] - noisefilter["R"][b]) < mindiff):
                                mindiff = abs(noisefilter["F"][a] - noisefilter["R"][b])
                            if (noisefilter["F"][a] < noisefilter["R"][b]):
                                a += 1
                            else:
                                b += 1
                        if k < mindiff <= (k+10):
                            cooccurencelist[runnum][k][n]["fr"] += 1
                            gaplist[runnum][k][n].append(mindiff-k)
                    elif seq1 in done and seq2 not in done:
                        cooccurencelist[runnum][k][n]["f"] += 1
                    elif seq2 in done and seq1 not in done:
                        cooccurencelist[runnum][k][n]["r"] += 1



def multiPosList(numofruns):
    for x in range(1, numofruns+1):
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
            avg5 = len(barcodeprimers53[runnum][0])
            avg3 = len(barcodeprimers53[runnum][1])
        except:
            continue
        c = 0
        if len(line) == l and "N" not in line:
            LSeqNums[runnum] += 1
            TSeqNums[runnum] += 1
            done = {}
            noisefilter = {"F": [],"R": []}
            for n in range(1, nmotifs+1):
                done[n] = set()
                for x in range(0,(l+1-k-avg5-avg3)):
                    kmers = kmer2hash(str(line[x+avg5:x+k+avg5]))
                    #if kmers not in done[n]:
                    if kmers < combinations:
                        if kmers in hamminglist2[runnum][k][n]:
                            poslist[runnum][k][n][x] += 1
                            done[n].add(kmers)
                        if kmers in rhamminglist2[runnum][k][n]:
                            rposlist[runnum][k][n][x] += 1
                            done[n].add(kmers)
                            if c == 0:
                                numoftfbsseq[runnum][k] += 1
                                c += 1
                    if x == int(l-k)-(avg5+avg3):
                        if seq1 & seq2 in done[n]:
                            numf = len(noisefilter["F"])
                            numr = len(noisefilter["R"])
                            a = 0
                            b = 0
                            mindiff = 1000
                            while (a < numf and b < numr):
                                if (k < abs(noisefilter["F"][a] - noisefilter["R"][b]) < mindiff):
                                    mindiff = abs(noisefilter["F"][a] - noisefilter["R"][b])
                                if (noisefilter["F"][a] < noisefilter["R"][b]):
                                    a += 1
                                else:
                                    b += 1
                            if k < mindiff < l:
                                cooccurencelist[runnum][k][n]["fr"] += 1
                                gaplist[runnum][k][n].append(mindiff-k)
                        elif seq1 in done and seq2 not in done[n]:
                            cooccurencelist[runnum][k][n]["f"] += 1
                        elif seq2 in done and seq1 not in done[n]:
                            cooccurencelist[runnum][k][n]["r"] += 1

        else:
            TSeqNums[runnum] += 1
            continue


def multiPosListmulti(numofruns):
    for j in files:
        fullname = os.path.join(datadir, j)
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
    for i in range(0, 6):
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
        bar.set_ylim(0,0.5)
        bar.set_yticks(yaxis)

        bar.set_xticks(xaxis)
        bar.bar(fxaxis, fseq, width, label = 'Forward strands')
        bar.bar(rxaxis, rseq, width, label = 'Reverse strands')
        bar.axhline(y=raverage, xmin=0.01, xmax=0.99, linestyle='dashed', color = 'black', linewidth = 0.75)

        bar.text((len(xaxis)-5), 0.4, "Average = "+str(round(raverage, 4)))
        handels=('Forward strands', 'Reverse strands')
        label=('Forward strands', 'Reverse strands')
        bar.legend(loc = 1)

    else:
        fig, kde = plt.subplots()
        kde.set_xlabel("TFBS relative position")
        kde.set_ylabel("Probability density function")
        kde.set_title("TFBS Position Distribution for Run: "+str(runnum+(startround-1))+", K: "+str(k)+", Consensus: "+str(consensuslist[runnum][k][n]))
        kde.set_xlim(0,1)
        x_d = np.linspace(0,1)

        flen = len(poslist[runnum][k][n])
        rlen = len(rposlist[runnum][k][n])
        tot = flen+rlen

        kdef = KernelDensity(bandwidth = 0.05, kernel='gaussian')
        kder = KernelDensity(bandwidth = 0.05, kernel='gaussian')

        kdef.fit(np.asarray(poslist[runnum][k][n])[:, None])
        kder.fit(np.asarray(rposlist[runnum][k][n])[:, None])

        logprobf = kdef.score_samples(x_d[:, None])
        logprobr = kder.score_samples(x_d[:, None])

        plt.fill_between(x_d, np.exp(logprobf*(flen/tot)), alpha=0.5)
        plt.fill_between(x_d, np.exp(logprobr*(rlen/tot)), alpha=0.5)

        label=("Forward", "Reverse")
        plt.legend(label)

    plt.savefig("figures/"+str(identifier)+"/position_bias/pos_"+str(identifier)+"_"+str(runnum+(startround-1))+"_"+str(k)+"_"+str(n), dpi=600)
    plt.close()



def plotrange(numofruns):
    for r in range(1, numofruns+1):
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
        for k in range(mink, maxk+1):
            numofkmers[i].update({k:sum(kmercount[i][k].values())})
            numofuniquekmers[i].update({k:len(kmercount[i][k])})
            numoftfbs[i].update({k:len(hamminglist2[i][k])})

numberofukmers()



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
                    seqbias[r][k].append(abs(seqbiasval))
                except:
                    continue
    return seqbias


seqbias = seqbiasfinder()


def cooccurenceindex():
    for r in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            for n in range(1, nmotifs+1):
                total = (cooccurencelist[r][k][n]["fr"]+cooccurencelist[r][k][n]["f"]+cooccurencelist[r][k][n]["r"])
                if total == 0:
                    total = 1
                cooccurenceindexlist[r][k][n] = (cooccurencelist[r][k][n]["fr"]/total)

cooccurenceindex()


def gapstdevcalc():
    for r in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            for n in range(1, nmotifs+1):
                if len(gaplist[r][k][n]) >= 2:
                    gapstdev[r][k][n]["mean"] += int(round((sum(gaplist[r][k][n])/len(gaplist[r][k][n])), 0))
                    gapstdev[r][k][n]["stdev"] += round(statistics.stdev(gaplist[r][k][n]),4)
                else:
                    gapstdev[r][k][n]["mean"] = 0
                    gapstdev[r][k][n]["stdev"] = 0

gapstdevcalc()


def writemetainfo():
    metafile = open("figures/"+str(identifier)+"/metainfo", "w")
    metafile.write("Analysis identifier: \t"+str(inputlist[0])+"\n")
    metafile.write("Data directory: \t"+str(inputlist[1])+"\n")
    metafile.write("Number of runs: \t"+str(inputlist[2])+"\n")
    metafile.write("Reverse compliment: \t"+str(inputlist[3])+"\n")
    metafile.write("Min k value: \t"+str(inputlist[4])+"\n")
    metafile.write("Max k value: \t"+str(inputlist[5])+"\n")
    metafile.write("Experiment type: \t"+str(inputlist[6])+"\n")
    metafile.write("Number of motifs: \t"+str(nmotifs)+"\n")
    metafile.write("Allowed hamming dist: \t"+str(allowham)+"\n")
    metafile.write("\n")
    if extype == "ChIP":
        for n in range(numofruns):
            metafile.write("File name "+str(n+1)+": \t"+str(inputlist[n+7])+"\n")
            metafile.write("Number of sequences "+str(n+1)+": \t"+str(TSeqNums[n+1])+"\n")
            if nobarcode == False:
                metafile.write("Barcode 5' "+str(n+1)+": \t"+str(barcodeprimers53[n+startround][0])+"\n")
                metafile.write("Barcode 3' "+str(n+1)+": \t"+str(barcodeprimers53[n+startround][1])+"\n")
            for k in range(mink, maxk+1):
                metafile.write("Sequence bias ("+str(n+1)+", "+str(k)+"): \t"+str(round((sum(seqbias[n+1][k])/len(seqbias[n+1][k])),7))+"\n")
                for x in range(1, nmotifs+1):
                    metafile.write("Consensus ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(consensuslist[n+1][k][x])+"\n")
                    metafile.write("PWM ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(logoform[n+1][k][x])+"\n")
                    metafile.write("Dimer co-occurence ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(cooccurenceindexlist[n+1][k][x])+"\n")
                    metafile.write("Mean gap ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(gapstdev[n+1][k][x]["mean"])+"\n")
                    metafile.write("Gap size stdev ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(gapstdev[n+1][k][x]["stdev"])+"\n")
            metafile.write("\n")
    if extype == "SELEX":
        for n in range(numofruns):
            metafile.write("File name "+str(n+1)+": \t"+str(inputlist[(n*2)+7])+"\n")
            metafile.write("Sequence length "+str(n+1)+": \t"+str(inputlist[(n*2)+8])+"\n")
            metafile.write("Total number of sequences "+str(n+1)+": \t"+str(TSeqNums[n+1])+"\n")
            metafile.write("Passed number of sequences "+str(n+1)+": \t"+str(LSeqNums[n+1])+"\n")
            if nobarcode == False:
                metafile.write("Barcode 5' "+str(n+1)+": \t"+str(barcodeprimers53[n+startround][0])+"\n")
                metafile.write("Barcode 3' "+str(n+1)+": \t"+str(barcodeprimers53[n+startround][1])+"\n")
            for k in range(mink, maxk+1):
                metafile.write("Sequence bias ("+str(n+1)+", "+str(k)+"): \t"+str(round((sum(seqbias[n+1][k])/len(seqbias[n+1][k])),7))+"\n")
                for x in range(1, nmotifs+1):
                    metafile.write("Consensus ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(consensuslist[n+1][k][x])+"\n")
                    metafile.write("PWM ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(logoform[n+1][k][x])+"\n")
                    metafile.write("Dimer co-occurence ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(cooccurenceindexlist[n+1][k][x])+"\n")
                    metafile.write("Mean gap ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(gapstdev[n+1][k][x]["mean"])+"\n")
                    metafile.write("Gap size stdev ("+str(n+1)+", "+str(k)+", "+str(x)+"): \t"+str(gapstdev[n+1][k][x]["stdev"])+"\n")
            metafile.write("\n")

    metafile.close()

writemetainfo()
