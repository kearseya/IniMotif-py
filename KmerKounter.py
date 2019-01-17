from itertools import groupby
from collections import Counter
from collections import OrderedDict
import numpy as np
#import sys
#import os


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
def seq2hash(kmer):
    k = len(kmer)
    base = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    kh = np.zeros(1,dtype='uint32')
    kh = base[kmer[0]]
    for tb in kmer[1:]:
        kh = kh<<2
        kh += base[tb]
    return kh
"""


revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev

identifier = str(input("Identifier:"))


numofruns = int(input("Number of runs:"))
#FileName = "NF1-2"
#l = 33
#runnum = 1
revcompwanted = bool(input("Reverse compliment:"))

"""
FileName = input("Fasta File Name:")
runnum = int(input("Run number:"))
l = int(input("Read lengths:"))
"""
#seqdict = ["NA",{},{},{},{},{},{},{},{},{},{}]
runlists = []
kmercount = []
hamlist = []
hamdict = []
pwm = []
filenames = {}
lvalues = {}

def dictinit(numofruns):
    for _ in range(numofruns + 1):
        runlists.append({})
        kmercount.append({})
        hamlist.append({})
        hamdict.append({})
        pwm.append({})


dictinit(numofruns)

"""
def sequencedictmaker(FileName, runnum, l):
    fastaFileName = open(FileName, "r")
    seqlist = []
    seqcount = {}
    for line in fastaFileName:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            barcode = line[1:]
            if barcode:
                continue
        seq = line
        if len(seq) == l:
            seqlist.append(seq2hash(seq))
    countseq = Counter(seqlist)
    lseqdict = {**countseq}
    return lseqdict

def addseqdict(FileName, runnum):
    seqdict[runnum].update(sequencedictmaker(FileName, runnum, l))

addseqdict(FileName, 1)
"""



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



def CreateKmerList(FileName, runnum, k):
    runlists[runnum][k] = []
    fastaFileName = open(FileName, "r")
    avg = barcodechecker(FileName)
    for line in fastaFileName:
        line = line.strip()
        if line.startswith(">"):
            continue
        if len(line) == l and "N" not in line:
            for x in range(0,((len(line)+1)-k)-(2*avg)):
                kmers = str(line[x+avg:x+k+avg])
                if len(kmers) > 0 and line.count(kmers) == 1:
                    runlists[runnum][k].append(kmer2hash(kmers))
        else:
            continue


        if revcompwanted == True:
            if len(line) == l and "N" not in line:
                for x in range(0,((len(line)+1)-k)-(2*avg)):
                    rkmers = revComp(line[x+avg:x+k+avg])
                    if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) == 1:
                        runlists[runnum][k].append(kmer2hash(rkmers))
            else:
                continue



def RangeKmerList(mink,maxk):
    for i in range(mink,maxk+1):
        CreateKmerList(FileName, runnum, i)

#RangeKmerList(6,8)



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

#KmerCounter()



def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer():
    for i in kmercount[runnum]:
        hamlist[runnum][i] = []
        #hconsensus = (list(kmercount[runnum][i].keys())[0])
        hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
        consensus = hash2kmer(hconsensus, i)
        for x in list(kmercount[runnum][i].keys()):
            values = hash2kmer(x,i)
            if hamming_distance(consensus, values) <= 1:
                hamlist[runnum][i].append(kmer2hash(values))

#listhammer()



def dicthammer():
    for i in hamlist[runnum]:
        hamdict[runnum][i] = {}
        test = { z : kmercount[runnum][i][z] for z in hamlist[runnum][i] }
        #M = sum(kmercount[runnum][i].values())
        M = sum(test.values())
        test = { z : (kmercount[runnum][i][z]/M) for z in test }
        hamdict[runnum][i].update(test)

#dicthammer()



def startpwm():
    for i in hamdict[runnum]:
        pwm[runnum][i] = {}
        for j in range(1,i+1):
            pwm[runnum][i][j] = {"A":0}
            pwm[runnum][i][j].update({"C":0})
            pwm[runnum][i][j].update({"G":0})
            pwm[runnum][i][j].update({"T":0})

#startpwm()



def pwmmaker():
    for i in hamdict[runnum]:
        for x in hamdict[runnum][i]:
            kmer = hash2kmer((x),i)
            for j in range(1,i+1):
                if kmer[j-1] == "A":
                    pwm[runnum][i][j]["A"] += hamdict[runnum][i][x]
                elif kmer[j-1] == "C":
                    pwm[runnum][i][j]["C"] += hamdict[runnum][i][x]
                elif kmer[j-1] == "T":
                    pwm[runnum][i][j]["T"] += hamdict[runnum][i][x]
                elif kmer[j-1] == "G":
                    pwm[runnum][i][j]["G"] += hamdict[runnum][i][x]

#pwmmaker()



def addrun():
    global FileName
    FileName = input("Fasta File Name:")
    global runnum
    runnum = int(input("Run number:"))
    filenames.update({runnum:FileName})
    global l
    l = int(input("Read lengths:"))
    lvalues.update({runnum:l})
    global mink
    mink = int(input("Minimum kmer:"))
    global maxk
    maxk = int(input("Maximum kmer:"))
    RangeKmerList(mink,maxk)
    KmerCounter()
    listhammer()
    dicthammer()
    startpwm()
    pwmmaker()



def addingall(n):
    for _ in range(n):
        addrun()


addingall(numofruns)


#print(runlists)
#print(kmercount)
#print(hamlist)
#print(hamdict)
#print(pwm)

print('Generating Logos')
import WebLogo
print('Generating Hamming distance figures')
import hamdistfig
print('Generating Position bias figures')
import postionbar
print('Generating Kmer frequency figures')
import kmerfreq
