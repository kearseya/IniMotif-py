import itertools
from collections import Counter
from collections import OrderedDict
import numpy as np
from Bio import SeqIO
#import sys
import os

#cwd = os.getcwd()

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
    kh = np.zeros(1,dtype='uint64')
    kh = base[kmer[0]]
    for tb in kmer[1:]:
        kh = kh<<2
        kh += base[tb]
    return kh
"""

bases = ['A', 'C', 'G', 'T']
revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

from GUI import inputlist
#print(inputlist)
from GUI import startround
from GUI import multiround
from GUI import knownbarcode

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev



def fiveprimerfinder(listin):
    primer = 0
    for i in range(len(listin[0])):
        truefor = 0
        for x in range(len(listin)):
            if listin[0][i:] == listin[x][i:]:
                truefor += 1
                if truefor == len(listin):
                    primer += 1
    return primer

def threeprimerfinder(listin):
    primer = 0
    for i in range(len(listin[0])):
        truefor = 0
        for x in range(len(listin)):
            if listin[0][i:] == listin[x][i:]:
                truefor += 1
                if truefor == len(listin):
                    primer += 1
    return primer



def initialinput():
    #print(inputlist)
    global identifier
    global datadir
    global numofruns
    global revcompwanted
    global mink
    global maxk
    global extype
    #try:
        #if len(inputlist) > 1:
    identifier = str(inputlist[0])
    datadir = inputlist[1]
    numofruns = int(inputlist[2])
    revcompwanted = inputlist[3]
    mink = int(inputlist[4])
    maxk = int(inputlist[5])
    extype = inputlist[6]

    if multiround == True:
        global files
        files = []
        for x in range(numofruns):
            if len(inputlist[(x*2)+7]) > 0:
                files.append(inputlist[(x*2)+7])
        #print("Files")
        #print(files)

    if knownbarcode == True:
        global barcodeprimers53
        barcodeprimers53 = {}
        for x in range(startround, (numofruns+startround)):
            barcodeprimers53[x] = (inputlist[(((x-startround)*2)+(numofruns*2)+7)], inputlist[(((x-startround)*2)+(numofruns*2)+8)])
        #print("barprim53")
        #print(barcodeprimers53)

        barcodeprimers5 = []
        for x in barcodeprimers53:
            barcodeprimers5.append(barcodeprimers53[x][0])
        #print("barcodeprimers5")
        #print(barcodeprimers5)
        len5 = []
        for i in barcodeprimers5:
            len5.append(len(i))
        #print("len5")
        #print(len5)
        max5 = max(len5)
        min5 = min(len5)
        if max5 == min5:
            primer5len = fiveprimerfinder(barcodeprimers5)
            #print("primer5len")
            #print(primer5len)
        barcodeprimers3 = []
        for x in barcodeprimers53:
            barcodeprimers3.append(barcodeprimers53[x][1])
        #print("barcodeprimers3")
        #print(barcodeprimers3)
        len3 = []
        for i in barcodeprimers3:
            len3.append(len(i))
        #print("len3")
        #print(len3)
        max3 = max(len3)
        min3 = min(len3)
        if max3 == min3:
            primer3len = threeprimerfinder(barcodeprimers3)
            #print("primer3len")
            #print(primer3len)

        global barcodes5
        barcodes5 = {}
        global barcodeslist5
        barcodeslist5 = []
        for x in barcodeprimers53:
            barcodes5[barcodeprimers53[x][0][:-primer5len]] = x
            barcodeslist5.append(barcodeprimers53[x][0][:-primer5len])
        #print("barcodes5")
        #print(barcodes5)
        #print("barcodeslist5")
        #print(barcodeslist5)
        list5 = []
        for i in barcodes5:
            list5.append(len(i))
        #print("list5")
        #print(list5)
        global barfiveslice
        barfiveslice = min(list5)
        #print("barfiveslice")
        #print(barfiveslice)
        global primer5
        primer5 = barcodeprimers5[0][primer5len:]
        #print("primer5")
        #print(primer5)
        global barcodes3
        barcodes3 = {}
        for x in barcodeprimers53:
            barcodes3[barcodeprimers53[x][1][-primer3len:]] = x
        #print("barcodes3")
        #print(barcodes3)
        list3 = []
        for i in barcodes3:
            list3.append(len(i))
        #print("list3")
        #print(list3)
        global barthreeslice
        barthreeslice = min(list3)
        #print("barthreeslice")
        #print(barthreeslice)
        global primer3
        primer3 = barcodeprimers3[0][:primer3len]
        #print("primer3")
        #print(primer3)

"""
    except:
        print("Command line input required")
        identifier = str(input("Identifier: "))
        datadir = input("Path to data directory: ")
        numofruns = int(input("Number of runs: "))
        revcompwanted = bool(input("Reverse compliment (any key, leave blank for False): "))
        extype = str(input("Experiment type: "))
        global startroundcl
        global multiroundcl
        global knownbarcodecl
        startroundcl = int(input("Start round: "))
        multiroundcl = bool(input("Multiple rounds per file(s): "))
        konwnbarcodecl = bool(input("Known barcodes/primers: "))
"""
initialinput()

"""
FileName = input("Fasta File Name:")
runnum = int(input("Run number:"))
l = int(input("Read lengths:"))
"""

#kmercombinations = {}
kmercount = []
hamlist = []
hamdict = []
pwm = []
filenames = {}
ufilenames = {}
filenames1 = {}
lvalues = {}

def dictinit(numofruns):
    for _ in range(numofruns + startround):
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
    while len(bar) <= 500:
        for line in fastaFileName:
            line = line.strip()
            if line.startswith(">") or line.startswith("@"):
                continue
            if len(line) == l and "N" not in line:
                for i in range(12):
                    start = line[0:i]
                    end = revComp(line[len(line)-i : len(line)])
                    if start == end:
                        x = i
                    else:
                        break
                bar.append(x)
            if line.startswith("+"):
                next(fastaFileName)
            else:
                continue
    favg = (sum(bar)/len(bar))
    avg = round(favg)
    return avg


"""
for hiding in hiding:
    def CreateKmerList(FileName, runnum, k):
        runlists[runnum][k] = []
        fastaFileName = open(FileName, "r")
        avg = barcodechecker(FileName)
        for line in fastaFileName:
            line = line.strip()
            if lineinsert.startswith(">"):
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
                        kmers = str(line[x+avg:x+k+avg])
                        rkmers = revComp(line[x+avg:x+k+avg])
                        if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) == 1:
                            runlists[runnum][k].append(kmer2hash(rkmers))
                else:
                    continue

    def RangeKmerList(mink,maxk):
        for i in range(mink,maxk+1):
            CreateKmerList(FileName, runnum, i)

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
"""


def Combinations(mink, maxk):
    for runnum in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            kmercount[runnum][k] = {}
            #kmercombinations[k] = [''.join(p)) for p in itertools.product(bases, repeat=k]
            for i in range(0, 4**k):
                kmercount[runnum][k][i] = 0

Combinations(mink, maxk)





def KmerCounterSELEX(FileName, runnum, k):
    fastaFileName = open(FileName, "r")
    if knownbarcode == False:
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
        if len(line) == l and "N" not in line:
            done = []
            for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                kmers = str(line[x+avg5:x+k+avg5])
                if kmers not in done:
                    if kmer2hash(kmers) < combinations:
                        kmercount[runnum][k][kmer2hash(kmers)] += 1
                        done.append(kmers)
                        #if revcompwanted == True:
                            #rkmers = revComp(kmers)
                            #if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                                #if kmer2hash(rkmers) < combinations:
                                    #kmercount[runnum][k][kmer2hash(rkmers)] += 1
        else:
            continue

        if revcompwanted == True:
            if len(line) == l and "N" not in line:
                rdone = []
                for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                    kmers = str(line[x+avg5:x+k+avg5])
                    rkmers = revComp(line[x+avg5:x+k+avg5])
                    if kmers not in rdone:
                        if kmer2hash(rkmers) < combinations:
                            kmercount[runnum][k][kmer2hash(rkmers)] += 1
                            rdone.append(rkmers)
            else:
                continue

def RangeKmerCounterSELEX(FileName, runnum, mink, maxk):
    for i in range(mink,maxk+1):
        KmerCounterSELEX(FileName, runnum, i)



def KmerCounterSELEXmulti(FileName, k):
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
        if len(line) == l and "N" not in line:
            if line[:barfiveslice] in barcodeslist5:
                done = []
                for x in range(0,(l+1-k-len(barcodeprimers53[runnum][0])-len(barcodeprimers53[runnum][1]))):
                    kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                    if kmers not in done:
                        if kmer2hash(kmers) < combinations:
                            kmercount[runnum][k][kmer2hash(kmers)] += 1
                            done.append(kmers)
                            #if revcompwanted == True:
                                #rkmers = revComp(kmers)
                                #if rkmers not in done:
                                    #if kmer2hash(rkmers) < combinations:
                                        #kmercount[runnum][k][kmer2hash(rkmers)] += 1
        else:
            continue

        if revcompwanted == True:
            if len(line) == l and "N" not in line:
                if line[:barfiveslice] in barcodeslist5:
                    runnum = barcodes5[line[:barfiveslice]]
                    rdone = []
                    for x in range(0,((len(line)+1)-k)-len(barcodeprimers53[runnum][0])-len(barcodeprimers53[runnum][1])):
                        kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                        rkmers = revComp(kmers)
                        if rkmers not in done:
                            if kmer2hash(rkmers) < combinations:
                                kmercount[runnum][k][kmer2hash(rkmers)] += 1
                                rdone.append(rkmers)
            else:
                continue

def RangeKmerCounterSELEXmulti(mink, maxk):
    for x in files:
        fullname = os.path.join(datadir, x)
        for i in range(mink,maxk+1):
            KmerCounterSELEXmulti(fullname, i)




def KmerCounterChipDNA(FileName, runnum, k):
    fastaFileName = open(FileName, "r")
    if knownbarcode == False:
        avg5 = barcodechecker(FileName)
        avg3 = avg5
    if knownbarcode == True:
        avg5 = len(barcodeprimers53[runnum][0])
        avg3 = len(barcodeprimers53[runnum][1])
    combinations = 4**k
    firstline = fastaFileName.readline()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
            kmers = str(line[x+avg5:x+k+avg5])
            try:
                if kmer2hash(kmers) < combinations:
                    kmercount[runnum][k][kmer2hash(kmers)] += 1
                    #if revcompwanted == True:
                        #rkmers = kmer2hash(revComp(kmers))
                        #kmercount[runnum][k][kmer2hash(rkmers)] += 1
            except:
                continue
        else:
            continue

        if revcompwanted == True:
            for x in range(0,((len(line)+1)-k-(avg5+avg3))):
                kmers = str(line[x+avg5:x+k+avg5])
                rkmers = revComp(line[x+avg5:x+k+avg5])
                if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                    if kmer2hash(rkmers) < combinations:
                        kmercount[runnum][k][kmer2hash(rkmers)] += 1
            else:
                continue

def RangeKmerCounterChipDNA(FileName, runnum, mink, maxk):
    for i in range(mink,maxk+1):
        KmerCounterChipDNA(FileName, runnum, i)


def KmerCounterChipDNAmulti(FileName, k):
    fastaFileName = open(FileName, "r")
    if knownbarcode == False:
        avg5 = barcodechecker(FileName)
        avg3 = avg5
    if knownbarcode == True:
        avg5 = len(barcodeprimers53[runnum][0])
        avg3 = len(barcodeprimers53[runnum][1])
    combinations = 4**k
    firstline = fastaFileName.readline()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        if line[:barfiveslice] in barcodeslist5:
            runnum = barcodes5[line[:barfiveslice]]
            for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                kmers = str(line[x+avg5:x+k+avg5])
                try:
                    if kmer2hash(kmers) < combinations:
                        kmercount[runnum][k][kmer2hash(kmers)] += 1
                        #if revcompwanted == True:
                            #rkmers = kmer2hash(revComp(kmers))
                            #kmercount[runnum][k][kmer2hash(rkmers)] += 1
                except:
                    continue
            else:
                continue

            if revcompwanted == True:
                for x in range(0,((len(line)+1)-k-(avg5+avg3))):
                    kmers = str(line[x+avg5:x+k+avg5])
                    rkmers = revComp(line[x+avg5:x+k+avg5])
                    if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                        if kmer2hash(rkmers) < combinations:
                            kmercount[runnum][k][kmer2hash(rkmers)] += 1
                else:
                    continue

def RangeKmerCounterChipDNAmulti(mink, maxk):
    for x in files:
        fullname = os.path.join(datadir, x)
        for i in range(mink,maxk+1):
            KmerCounterChipDNAmulti(fullname, i)





def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))



def listhammer(runnum):
    for i in kmercount[runnum]:
        hamlist[runnum][i] = []
        #hconsensus = (list(kmercount[runnum][i].keys())[0])
        hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
        consensus = hash2kmer(hconsensus, i)
        for x in list(kmercount[runnum][i].keys()):
            values = hash2kmer(x,i)
            if hamming_distance(consensus, values) <= 1:
                hamlist[runnum][i].append(x)



def dicthammer(runnum):
    for i in hamlist[runnum]:
        hamdict[runnum][i] = {}
        test = { z : kmercount[runnum][i][z] for z in hamlist[runnum][i] }
        #M = sum(kmercount[runnum][i].values())
        M = sum(test.values())
        test = { z : (kmercount[runnum][i][z]/M) for z in test }
        hamdict[runnum][i].update(test)




def startpwm(runnum):
    for i in hamdict[runnum]:
        pwm[runnum][i] = {}
        for j in range(1,i+1):
            pwm[runnum][i][j] = {"A":0}
            pwm[runnum][i][j].update({"C":0})
            pwm[runnum][i][j].update({"G":0})
            pwm[runnum][i][j].update({"T":0})



def pwmmaker(runnum):
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



def addruncl(runnum):
    global mutliround
    multiround = bool(input("Multiple roudns per file: "))
    #if multiround == True:
    #    global numoffiles
    #    numoffiles = int(input("Number of files: "))
    #    for _ in range(numoffiles):
    #        global FileName1
    #        global FileName
    #        FileName1 = input("Fasta File Name:")
    #        FileName = os.path.join(datadir, FileName1)


    if multiround == False:
        global FileName1
        global FileName
        FileName1 = input("Fasta File Name:")
        FileName = os.path.join(datadir, FileName1)
        #global runnum
        #runnum = int(input("Run number:"))
        filenames1.update({runnum:FileName1})
        ufilenames.update({FileName:runnum})
        filenames.update({runnum:FileName})
        global l
        l = int(input("Read lengths:"))
        lvalues.update({runnum:l})
        #global mink
        #mink = int(input("Minimum kmer:"))
        #global maxk
        #maxk = int(input("Maximum kmer:"))
        #Combinations(runnum, mink, maxk)
        RangeKmerCounter(FileName, runnum, mink, maxk)
        listhammer(runnum)
        dicthammer(runnum)
        startpwm(runnum)
        pwmmaker(runnum)

def addrungui():
    try:
        if len(inputlist) > 1:
            global FileName1
            global FileName
            global runnum
            global l
            if multiround == False:
                if extype == "SELEX":
                    for x in range(0, numofruns):
                        FileName1 = str(inputlist[(x*2)+7])
                        FileName = os.path.join(datadir, FileName1)
                        runnum = x+startround
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        l = int(inputlist[(x*2)+8])
                        lvalues.update({runnum:l})
                        RangeKmerCounterSELEX(FileName, runnum, mink, maxk)
                        listhammer(runnum)
                        dicthammer(runnum)
                        startpwm(runnum)
                        pwmmaker(runnum)
                if extype in ["ChIP", "DNase"]:
                    for x in range(0, numofruns):
                        FileName1 = str(inputlist[x+7])
                        FileName = os.path.join(datadir, FileName1)
                        runnum = x+1
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        #Combinations(runnum, mink, maxk)
                        RangeKmerCounterChipDNA(FileName, runnum, mink, maxk)
                        listhammer(runnum)
                        dicthammer(runnum)
                        startpwm(runnum)
                        pwmmaker(runnum)
            if multiround == True:
                if extype == "SELEX":
                    for x in range(startround, (numofruns+startround)):
                        l = int(inputlist[((x-startround)*2)+8])
                        lvalues.update({x:l})
                    #l = int(inputlist[8])
                    RangeKmerCounterSELEXmulti(mink, maxk)
                    for x in range(startround, (numofruns+startround)):
                        listhammer(x)
                        dicthammer(x)
                        startpwm(x)
                        pwmmaker(x)
                if extype in ["ChIP", "DNase"]:
                    RangeKmerCounterChipDNAmulti(mink, maxk)
                    for x in range(startround, numofrounds):
                        listhammer(x)
                        dicthammer(x)
                        startpwm(x)
                        pwmmaker(x)


    except:
        print("Command line input required (removed)")


def addingall(n):
    try:
        if len(inputlist) > 1:
            addrungui()
    except:
        for x in range(1, n+1):
            #addruncl(x)
            print("Removed")


addingall(numofruns)


#print(inputlist)


def removezeros():
    for x in range(1, numofruns+1):
         for k in range(mink, maxk+1):
             for i in range(0, 4**k):
                 if kmercount[x][k][i] == 0:
                     del kmercount[x][k][i]

removezeros()



#print(kmercombinations)

#print("kmercount")
#print(kmercount)
#print("hamlist")
#print(hamlist)
#print("hamdict")
#print(hamdict)
#print("pwm")
#print(pwm)


print('Generating Logos')
import WebLogoMod

print('Generating Hamming distance figures')
import HamDistFig

if extype == "SELEX" or "ATAC":
    print('Generating Position bias figures')
    import PositionBias

if numofruns > 1:
    print('Generating Kmer frequency figures')
    import KmerFreq
