import itertools
from collections import Counter
from collections import OrderedDict
import numpy as np
from Bio import SeqIO
#import sys
import os
import readline

import operator

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

def importfromgui():
    global startround
    global multiround
    global knownbarcode
    global nobarcode
    global inputlist
    try:
        from GUI import inputlist
        #print(inputlist)
        from GUI import startround
        from GUI import multiround
        from GUI import knownbarcode
        from GUI import nobarcode
    except:
        print("Command line input required")
        inputlist = ["cl"]
        multiround = str(input("Multiple rounds per file (y/n)?: "))
        knownbarcode = str(input("Known barcodes (y/n)?: "))
        startround = int(input("Starting round: "))
        if multiround in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
            multiround = True
        else:
            multiround = False
        if knownbarcode in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
            nobarcode = False
            knownbarcode = True
        elif knownbarcode in ["NA", "na", "None", "none"]:
            nobarcode = True
            knownbarcode = False
        else:
            nobarcode = False
            knownbarcode = False

importfromgui()


def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev


"""
def fiveprimerfinder(listin):
    primer = 0
    for i in range(len(listin[0])):
        truefor = 0
        #print(diff)
        for x in range(len(listin)):
            diff = abs((len(listin[0])-len(listin[x])))
            if diff == 0:
                if listin[0][i:] == listin[x][i:]:
                    truefor += 1
                    if truefor == len(listin):
                        primer += 1
            if diff != 0:
                if listin[0][i:-diff] == listin[x][i:-diff]:
                    truefor += 1
                    if truefor == len(listin):
                        primer += 1
    return primer
"""
def fiveprimerfinder(listin, min):
    primer = 0
    if nobarcode == True:
        return primer
    for i in range(min+1):
        truefor = 0
        for x in range(len(listin)):
            if listin[0][i:] == listin[x][i:]:
                truefor += 1
                if truefor == len(listin):
                    primer += 1
    return primer

def threeprimerfinder(listin, min):
    primer = 0
    if nobarcode == True:
        return primer
    for i in range(min+1):
        truefor = 0
        for x in range(len(listin)):
            if listin[0][i:] == listin[x][i:]:
                truefor += 1
                if truefor == len(listin):
                    primer += 1
    return primer


filenames = {}
ufilenames = {}
filenames1 = {}
lvalues = {}

def input_with_prefill(prompt, text):
    def hook():
        readline.insert_text(text)
        readline.redisplay()
    readline.set_pre_input_hook(hook)
    result = input(prompt)
    readline.set_pre_input_hook()
    return result

def autolenfinder(run):
    cwf = open(filenames[run])
    for linenum, line in enumerate(cwf):
        if linenum % 4 == 1:
            if "N" not in line:
                firstlval = len(line.strip())
                return firstlval
        if linenum > 32:
            break

def autofiveprimefinder(run, fivesplice):
    cwf = open(filenames[run])
    if nobarcode == True:
        five = 0
        return five
    for linenum, line in enumerate(cwf):
        if linenum % 4 == 1:
            line = line.strip()
            if extype == "SELEX":
                if len(line) == int(lvalues[run]) and "N" not in line:
                    five = line[:fivesplice]
                    print("5': "+str(five))
                    return five
            else:
                if "N" not in line:
                    five = line[:fivesplice]
                    print("5': NA "+str(five))
                    return five
        if linenum > 32:
            break

def autothreeprimefinder(run, threesplice):
    cwf = open(filenames[run])
    if nobarcode == True:
        three = 0
        return three
    for linenum, line in enumerate(cwf):
        if linenum % 4 == 1:
            line = line.strip()
            if extype == "SELEX":
                if len(line) == int(lvalues[run]) and "N" not in line:
                    three = line[-threesplice:]
                    print("3': "+str(three))
                    return three
            else:
                if "N" not in line:
                    three = line[-threesplice:]
                    print("3': "+str(three))
                    return three
        if linenum > 32:
            break


def barcodechecker(FileName):
    global nobarcode
    if nobarcode == True:
        avg = 0
        return avg
    bar = []
    fastaFileName = open(FileName, "r")
    while len(bar) <= 500:
        for line in fastaFileName:
            line = line.strip()
            if line.startswith(">") or line.startswith("@"):
                continue
            if extype == "SELEX":
                if len(line) == l and "N" not in line:
                    for i in range(4, 12):
                        start = line[0:i]
                        end = revComp(line[len(line)-i : len(line)])
                        if start == end:
                            bar.append(i)
                        else:
                            continue
                    if not bar:
                        avg = 0
                        nobarcode = True
                        return avg
                if line.startswith("+"):
                    next(fastaFileName)
                else:
                    continue
            else:
                if "N" not in line:
                    for i in range(4, 12):
                        start = line[0:i]
                        end = revComp(line[len(line)-i : len(line)])
                        if start == end:
                            bar.append(i)
                        else:
                            continue
                    if not bar:
                        avg = 0
                        nobarcode = True
                        return avg
                if line.startswith("+"):
                    next(fastaFileName)
                else:
                    continue
    favg = (sum(bar)/len(bar))
    avg = round(favg)
    if not avg:
        nobarcode = True
    return avg


def initialinput():
    #print(inputlist)
    global identifier
    global datadir
    global numofruns
    global revcompwanted
    global mink
    global maxk
    global extype
    #global multiround
    #global knownbarcode
    global barcodeprimers53
    barcodeprimers53 = {}
    global barcodes5
    global barcodeslist5
    global barfiveslice
    global primer5
    global barcodes3
    global barthreeslice
    global primer3

    if len(inputlist) > 1:
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
            #if max5 == min5:
            primer5len = fiveprimerfinder(barcodeprimers5, min5)
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
            #if max3 == min3:
            primer3len = threeprimerfinder(barcodeprimers3, min3)
            #print("primer3len")
            #print(primer3len)
            barcodes5 = {}
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
            barfiveslice = min(list5)
            #print("barfiveslice")
            #print(barfiveslice)
            primer5 = barcodeprimers5[0][primer5len:]
            #print("primer5")
            #print(primer5)
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
            barthreeslice = min(list3)
            #print("barthreeslice")
            #print(barthreeslice)
            primer3 = barcodeprimers3[0][:primer3len]
            #print("primer3")
            #print(primer3)

        print("Counting kmers")

    else:
        identifier = str(input("Identifier: "))
        datadir = input("Path to data directory: ")
        numofruns = int(input("Number of runs: "))
        revcompwanted = str(input("Reverse compliment (y/n)?: "))
        if revcompwanted in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
            revcompwanted = True
        else:
            revcompwanted = False
        extype = str(input("Experiment type (SELEX, ChIP, ATAC): "))
        mink = int(input("Minimum kmer length: "))
        maxk = int(input("Maximum kmer length: "))
        global clnmotifs
        global cllogotype
        global clallowham
        clnmotifs = int(input("Number of motifs: "))
        cllogotype = input("Bits (b) or Frequency (f)?: ")
        clallowham = int(input("Allowed hamming distance: "))

        global tryauto

        try:
            namesindirectory = os.listdir(str(datadir))
            orderednumbers = [99999999]
            for n in namesindirectory:
                if n[:3].isalpha() and n[3:-6].isnumeric():
                    orderednumbers.append(int(n[3:-6]))
            orderednumbers.sort(reverse=True)
            sixnumbersstart = str(orderednumbers[startround])
            if len(sixnumbersstart) < 6:
                sixnumbersstart = (6-len(sixnumbersstart))*"0"+sixnumbersstart
            for f in namesindirectory:
                if sixnumbersstart in f:
                    startingfile = f
        except:
            tryauto = False


        global FileName1
        global FileName



        if multiround == False:
            if extype == "SELEX":
                global l
                FileName1 = str(input("Fastx File Name (try \"auto\"): "))
                if FileName1 == "auto":
                    tryauto = True
                    FileName1 = str(input_with_prefill("Fastx File Name: ", startingfile))
                    FileName = os.path.join(datadir, FileName1)
                    #global runnum
                    #runnum = int(input("Run number:"))
                    filenames1.update({1:FileName1})
                    ufilenames.update({FileName:1})
                    filenames.update({1:FileName})
                    l = int(input_with_prefill("Read lengths: ", str(autolenfinder(1))))
                    lvalues.update({1:l})
                else:
                    FileName = os.path.join(datadir, FileName1)
                    #global runnum
                    #runnum = int(input("Run number:"))
                    filenames1.update({1:FileName1})
                    ufilenames.update({FileName:1})
                    filenames.update({1:FileName})
                    l = int(input("Read lengths: "))
                    lvalues.update({1:l})
                if knownbarcode == True:
                    fivein = str(input("5' barcode/primer ("+str(1+startround-1)+"): "))
                    if multiround == False and fivein.isnumeric():
                        fivein = autofiveprimefinder(1, int(fivein))
                    threein = str(input("3' primer/barcode ("+str(1+startround-1)+"): "))
                    if multiround == False and threein.isnumeric():
                        threein = autothreeprimefinder(1, int(threein))
                    barcodeprimers53[1] = (fivein, threein)
                if knownbarcode == False:
                    bothslice = barcodechecker(filenames[1])
                    fivein = autofiveprimefinder(1, int(bothslice))
                    threein = autothreeprimefinder(1, int(bothslice))
                    barcodeprimers53[1] = (fivein, threein)
                if tryauto == True:
                    for runnum in range(2, numofruns+1):
                        for f in namesindirectory:
                            if str(orderednumbers[startround+runnum-1]) in f:
                                nextfile = f
                        FileName1 = str(input_with_prefill("Fastx File Name: ", nextfile))
                        FileName = os.path.join(datadir, FileName1)
                        #global runnum
                        #runnum = int(input("Run number:"))
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        #global l
                        l = int(input_with_prefill("Read lengths: ", str(autolenfinder(2))))
                        lvalues.update({runnum:l})
                        if knownbarcode == True:
                            fivein = str(input("5' barcode/primer ("+str(runnum+startround-1)+"): "))
                            if multiround == False and fivein.isnumeric():
                                fivein = autofiveprimefinder(runnum, int(fivein))
                            threein = str(input("3' primer/barcode ("+str(runnum+startround-1)+"): "))
                            if multiround == False and threein.isnumeric():
                                threein = autothreeprimefinder(runnum, int(threein))
                            barcodeprimers53[runnum] = (fivein, threein)
                        if knownbarcode == False:
                            bothslice = barcodechecker(filenames[1])
                            fivein = autofiveprimefinder(1, int(bothslice))
                            threein = autothreeprimefinder(1, int(bothslice))
                            barcodeprimers53[1] = (fivein, threein)
                else:
                    for runnum in range(2, numofruns+1):
                        FileName1 = str(input("Fastx File Name: "))
                        FileName = os.path.join(datadir, FileName1)
                        #global runnum
                        #runnum = int(input("Run number:"))
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        #global l
                        l = int(input("Read lengths: "))
                        lvalues.update({runnum:l})
                        if knownbarcode == True:
                            fivein = str(input("5' barcode/primer ("+str(runnum+startround-1)+"): "))
                            if multiround == False and fivein.isnumeric():
                                fivein = autofiveprimefinder(runnum, int(fivein))
                            threein = str(input("3' primer/barcode ("+str(runnum+startround-1)+"): "))
                            if multiround == False and threein.isnumeric():
                                threein = autothreeprimefinder(runnum, int(threein))
                            barcodeprimers53[runnum] = (fivein, threein)
                        if knownbarcode == False:
                            bothslice = barcodechecker(filenames[runnum])
                            fivein = autofiveprimefinder(runnum, int(bothslice))
                            threein = autothreeprimefinder(runnum, int(bothslice))
                            barcodeprimers53[runnum] = (fivein, threein)
            else:
                FileName1 = str(input("Fastx File Name (try \"auto\"): "))
                if FileName1 == "auto":
                    tryauto = True
                    FileName1 = str(input_with_prefill("Fastx File Name: ", startingfile))
                    FileName = os.path.join(datadir, FileName1)
                    #global runnum
                    #runnum = int(input("Run number:"))
                    filenames1.update({1:FileName1})
                    ufilenames.update({FileName:1})
                    filenames.update({1:FileName})
                else:
                    FileName = os.path.join(datadir, FileName1)
                    #global runnum
                    #runnum = int(input("Run number:"))
                    filenames1.update({1:FileName1})
                    ufilenames.update({FileName:1})
                    filenames.update({1:FileName})
                if knownbarcode == True:
                    fivein = str(input("5' barcode/primer ("+str(1+startround-1)+"): "))
                    if multiround == False and fivein.isnumeric():
                        fivein = autofiveprimefinder(1, int(fivein))
                    threein = str(input("3' primer/barcode ("+str(1+startround-1)+"): "))
                    if multiround == False and threein.isnumeric():
                        threein = autothreeprimefinder(1, int(threein))
                    barcodeprimers53[1] = (fivein, threein)
                if knownbarcode == False:
                    bothslice = barcodechecker(filenames[1])
                    fivein = autofiveprimefinder(1, int(bothslice))
                    threein = autothreeprimefinder(1, int(bothslice))
                    barcodeprimers53[1] = (fivein, threein)
                if tryauto == True:
                    for runnum in range(2, numofruns+1):
                        for f in namesindirectory:
                            if str(orderednumbers[startround+runnum-1]) in f:
                                nextfile = f
                        FileName1 = str(input_with_prefill("Fastx File Name: ", nextfile))
                        FileName = os.path.join(datadir, FileName1)
                        #global runnum
                        #runnum = int(input("Run number:"))
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        if knownbarcode == True:
                            fivein = str(input("5' barcode/primer ("+str(runnum+startround-1)+"): "))
                            if multiround == False and fivein.isnumeric():
                                fivein = autofiveprimefinder(runnum, int(fivein))
                            threein = str(input("3' primer/barcode ("+str(runnum+startround-1)+"): "))
                            if multiround == False and threein.isnumeric():
                                threein = autothreeprimefinder(runnum, int(threein))
                            barcodeprimers53[runnum] = (fivein, threein)
                        if knownbarcode == False:
                            bothslice = barcodechecker(filenames[runnum])
                            fivein = autofiveprimefinder(runnum, int(bothslice))
                            threein = autothreeprimefinder(runnum, int(bothslice))
                            barcodeprimers53[runnum] = (fivein, threein)
                else:
                    for runnum in range(2, numofruns+1):
                        FileName1 = str(input("Fastx File Name: "))
                        FileName = os.path.join(datadir, FileName1)
                        #global runnum
                        #runnum = int(input("Run number:"))
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        if knownbarcode == True:
                            fivein = str(input("5' barcode/primer ("+str(runnum+startround-1)+"): "))
                            if multiround == False and fivein.isnumeric():
                                fivein = autofiveprimefinder(runnum, int(fivein))
                            threein = str(input("3' primer/barcode ("+str(runnum+startround-1)+"): "))
                            if multiround == False and threein.isnumeric():
                                threein = autothreeprimefinder(runnum, int(threein))
                            barcodeprimers53[runnum] = (fivein, threein)
                        if knownbarcode == False:
                            bothslice = barcodechecker(filenames[runnum])
                            fivein = autofiveprimefinder(runnum, int(bothslice))
                            threein = autothreeprimefinder(runnum, int(bothslice))
                            barcodeprimers53[runnum] = (fivein, threein)

        if multiround == True:
            global numoffiles
            numoffiles = int(input("Number of files: "))
            global clfiles
            clfiles = []
            for _ in range(numoffiles):
                FileName1 = input("Fastx File Name: ")
                FileName = os.path.join(datadir, FileName1)
                clfiles.append(FileName)
            files = clfiles
            if extype == "SELEX":
                for x in range(0, numofruns):
                    l = int(input("Length of reads ("+str(x+1)+"): "))
                    lvalues.update({(x+1):l})
                    fivein = str(input("5' barcode/primer ("+str(x+startround)+"): "))
                    threein = str(input("3' primer/barcode ("+str(x+startround)+"): "))
                    barcodeprimers53[x+1] = (fivein, threein)
            else:
                for x in range(1, numofruns+1):
                    fivein = str(input("5' barcode/primer ("+str(x+startround-1)+"): "))
                    threein = str(input("3' primer/barcode ("+str(x+startround-1)+"): "))
                    barcodeprimers53[x] = (fivein, threein)


        print("Counting kmers")

        if knownbarcode == True:
            """
            barcodeprimers53 = {}
            for x in range(startround, (numofruns+startround)):
                fivein = str(input("5' barcode/primer ("+str(x)+"): "))
                threein = str(input("3' primer/barcode ("+str(x)+"): "))
                barcodeprimers53[x] = (fivein, threein)
            """
            barcodeprimers5 = []
            for x in barcodeprimers53:
                barcodeprimers5.append(barcodeprimers53[x][0])
            len5 = []
            for i in barcodeprimers5:
                len5.append(len(i))
            max5 = max(len5)
            min5 = min(len5)
            #if max5 == min5:
            primer5len = fiveprimerfinder(barcodeprimers5, min5)
            barcodeprimers3 = []
            for x in barcodeprimers53:
                barcodeprimers3.append(barcodeprimers53[x][1])
            len3 = []
            for i in barcodeprimers3:
                len3.append(len(i))
            max3 = max(len3)
            min3 = min(len3)
            #if max3 == min3:
            primer3len = threeprimerfinder(barcodeprimers3, min3)
            barcodes5 = {}
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
            barfiveslice = min(list5)
            primer5 = barcodeprimers5[0][primer5len:]
            barcodes3 = {}
            for x in barcodeprimers53:
                barcodes3[barcodeprimers53[x][1][-primer3len:]] = x
            list3 = []
            for i in barcodes3:
                list3.append(len(i))
            barthreeslice = min(list3)
            primer3 = barcodeprimers3[0][:primer3len]


initialinput()

"""
FileName = input("Fastx File Name:")
runnum = int(input("Run number:"))
l = int(input("Read lengths:"))
"""

#kmercombinations = {}
kmercount = []
#hamlist = []
#hamdict = []
#pwm = []


def dictinit(numofruns):
    for _ in range(numofruns + startround):
        kmercount.append({})
        #hamlist.append({})
        #hamdict.append({})
        #pwm.append({})


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
    l = lvalues[runnum]
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
            done = set()
            for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                kmers = str(line[x+avg5:x+k+avg5])
                if kmers not in done:
                    if kmer2hash(kmers) < combinations:
                        kmercount[runnum][k][kmer2hash(kmers)] += 1
                        done.add(kmers)
                        #if revcompwanted == True:
                            #rkmers = revComp(kmers)
                            #if len(rkmers) > 0 and (line.count(kmers) + line.count(rkmers)) <= 2:
                                #if kmer2hash(rkmers) < combinations:
                                    #kmercount[runnum][k][kmer2hash(rkmers)] += 1
        else:
            continue

        if revcompwanted == True:
            if len(line) == l and "N" not in line:
                rdone = set()
                for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
                    kmers = str(line[x+avg5:x+k+avg5])
                    rkmers = revComp(line[x+avg5:x+k+avg5])
                    if kmers not in rdone:
                        if kmer2hash(rkmers) < combinations:
                            kmercount[runnum][k][kmer2hash(rkmers)] += 1
                            rdone.add(rkmers)
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
                done = set()
                for x in range(0,(l+1-k-len(barcodeprimers53[runnum][0])-len(barcodeprimers53[runnum][1]))):
                    kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                    if kmers not in done:
                        if kmer2hash(kmers) < combinations:
                            kmercount[runnum][k][kmer2hash(kmers)] += 1
                            done.add(kmers)
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
                    rdone = set()
                    for x in range(0,((len(line)+1)-k)-len(barcodeprimers53[runnum][0])-len(barcodeprimers53[runnum][1])):
                        kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                        rkmers = revComp(kmers)
                        if rkmers not in done:
                            if kmer2hash(rkmers) < combinations:
                                kmercount[runnum][k][kmer2hash(rkmers)] += 1
                                rdone.add(rkmers)
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
        done = set()
        for x in range(0,((len(line)+1)-k)-(avg5+avg3)):
            kmers = str(line[x+avg5:x+k+avg5])
            try:
                if kmers not in done:
                    if kmer2hash(kmers) < combinations:
                        kmercount[runnum][k][kmer2hash(kmers)] += 1
                        done.add(kmers)
                        #if revcompwanted == True:
                            #rkmers = kmer2hash(revComp(kmers))
                            #kmercount[runnum][k][kmer2hash(rkmers)] += 1
            except:
                continue
        else:
            continue

        if revcompwanted == True:
            rdone = set()
            for x in range(0,((len(line)+1)-k-(avg5+avg3))):
                kmers = str(line[x+avg5:x+k+avg5])
                rkmers = revComp(line[x+avg5:x+k+avg5])
                if len(rkmers) > 0 and rkmers not in rdone:
                    if kmer2hash(rkmers) < combinations:
                        kmercount[runnum][k][kmer2hash(rkmers)] += 1
                        rdone.add(rkmers)
            else:
                continue

def RangeKmerCounterChipDNA(FileName, runnum, mink, maxk):
    for i in range(mink,maxk+1):
        KmerCounterChipDNA(FileName, runnum, i)


def KmerCounterChipDNAmulti(FileName, k):
    fastaFileName = open(FileName, "r")
    combinations = 4**k
    firstline = fastaFileName.readline()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    for sequence in SeqIO.parse(FileName, filetype):
        line = str(sequence.seq)
        if line[:barfiveslice] in barcodeslist5:
            done = set()
            runnum = barcodes5[line[:barfiveslice]]
            for x in range(0,((len(line)+1)-k)-(len(barcodeprimers53[runnum][0])+len(barcodeprimers53[runnum][1]))):
                kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                try:
                    if kmers not in done:
                        if kmer2hash(kmers) < combinations:
                            kmercount[runnum][k][kmer2hash(kmers)] += 1
                            done.add(kmers)
                            #if revcompwanted == True:
                                #rkmers = kmer2hash(revComp(kmers))
                                #kmercount[runnum][k][kmer2hash(rkmers)] += 1
                except:
                    continue
            else:
                continue

            if revcompwanted == True:
                rdone = set()
                for x in range(0,((len(line)+1)-k-(len(barcodeprimers53[runnum][0])+len(barcodeprimers53[runnum][1])))):
                    kmers = str(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                    rkmers = revComp(line[x+len(barcodeprimers53[runnum][0]):x+k+len(barcodeprimers53[runnum][0])])
                    if len(rkmers) > 0 and rkmers not in rdone:
                        if kmer2hash(rkmers) < combinations:
                            kmercount[runnum][k][kmer2hash(rkmers)] += 1
                            rdone.add(rkmers)
                else:
                    continue

def RangeKmerCounterChipDNAmulti(mink, maxk):
    for x in files:
        fullname = os.path.join(datadir, x)
        for i in range(mink,maxk+1):
            KmerCounterChipDNAmulti(fullname, i)




"""
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


"""


def addruncl():
    #global FileName1
    #global FileName
    if multiround == False:
        if extype == "SELEX":
            """
            for runnum in range(1, numofruns+1):
                FileName1 = str(input("Fastx File Name: "))
                FileName = os.path.join(datadir, FileName1)
                #global runnum
                #runnum = int(input("Run number:"))
                filenames1.update({runnum:FileName1})
                ufilenames.update({FileName:runnum})
                filenames.update({runnum:FileName})
                global l
                l = int(input("Read lengths: "))
                lvalues.update({runnum:l})
            """
            for run in range(1, numofruns+1):
                RangeKmerCounterSELEX(filenames[run], run, mink, maxk)
                #listhammer(runnum)
                #dicthammer(runnum)
                #startpwm(runnum)
                #pwmmaker(runnum)
        if extype in ["ChIP", "DNase"]:
            for runnum in range(1, numofruns+1):
                RangeKmerCounterChipDNA(FileName, runnum, mink, maxk)
                #listhammer(runnum)
                #dicthammer(runnum)
                #startpwm(runnum)
                #pwmmaker(runnum)
    if multiround == True:
        """
        global numoffiles
        numoffiles = int(input("Number of files: "))
        global files
        files = []
        for _ in range(numoffiles):
            FileName1 = input("Fastx File Name: ")
            FileName = os.path.join(datadir, FileName1)
            files.append(FileName)
        """
        if extype == "SELEX":
            """
            for x in range(0, numofruns):
                l = int(input("Length of reads ("+str(x+1)+"): "))
                lvalues.update({(x+1):l})
            """
            RangeKmerCounterSELEXmulti(mink, maxk)
            #for x in range(1, numofruns+1):
                #listhammer(x)
                #dicthammer(x)
                #startpwm(x)
                #pwmmaker(x)
        if extype in ["ChIP", "DNase"]:
            RangeKmerCounterChipDNAmulti(mink, maxk)
            #for x in range(1, numofrounds+1):
                #listhammer(x)
                #dicthammer(x)
                #startpwm(x)
                #pwmmaker(x)

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
                        runnum = x+1
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        l = int(inputlist[(x*2)+8])
                        lvalues.update({runnum:l})
                        if knownbarcode == False:
                            bothslice = barcodechecker(filenames[runnum])
                            fivein = autofiveprimefinder(runnum, int(bothslice))
                            threein = autothreeprimefinder(runnum, int(bothslice))
                            barcodeprimers53[runnum] = (fivein, threein)
                        RangeKmerCounterSELEX(FileName, runnum, mink, maxk)
                        #listhammer(runnum)
                        #dicthammer(runnum)
                        #startpwm(runnum)
                        #pwmmaker(runnum)
                if extype in ["ChIP", "DNase"]:
                    for x in range(0, numofruns):
                        FileName1 = str(inputlist[x+7])
                        FileName = os.path.join(datadir, FileName1)
                        runnum = x+1
                        filenames1.update({runnum:FileName1})
                        ufilenames.update({FileName:runnum})
                        filenames.update({runnum:FileName})
                        if knownbarcode == False:
                            bothslice = barcodechecker(filenames[runnum])
                            fivein = autofiveprimefinder(runnum, int(bothslice))
                            threein = autothreeprimefinder(runnum, int(bothslice))
                            barcodeprimers53[runnum] = (fivein, threein)
                        #Combinations(runnum, mink, maxk)
                        RangeKmerCounterChipDNA(FileName, runnum, mink, maxk)
                        #listhammer(runnum)
                        #dicthammer(runnum)
                        #startpwm(runnum)
                        #pwmmaker(runnum)
            if multiround == True:
                if extype == "SELEX":
                    for x in range(0, numofruns):
                        l = int(inputlist[(x*2)+8])
                        lvalues.update({(x+1):l})
                    #l = int(inputlist[8])
                    RangeKmerCounterSELEXmulti(mink, maxk)
                    #for x in range(1, numofruns+1):
                        #listhammer(x)
                        #dicthammer(x)
                        #startpwm(x)
                        #pwmmaker(x)
                if extype in ["ChIP", "DNase"]:
                    RangeKmerCounterChipDNAmulti(mink, maxk)
                    #for x in range(1, numofrounds+1):
                        #listhammer(x)
                        #dicthammer(x)
                        #startpwm(x)
                        #pwmmaker(x)


    except:
        print("Further command line input required")


def addingall(n):
    if len(inputlist) > 1:
        print("adding")
        addrungui()
    else:
        addruncl()


addingall(numofruns)

"""
print("kmercount")
print(kmercount)
print("hamlist")
print(hamlist)
print("hamdict")
print(hamdict)
print("pwm")
print(pwm)
"""
#print(inputlist)


def removezeros():
    for x in range(1, numofruns+1):
         for k in range(mink, maxk+1):
             for i in range(0, 4**k):
                 if kmercount[x][k][i] == 0:
                     del kmercount[x][k][i]

removezeros()

top12t = []
top12 = []

def top12maker(numofruns):
    for _ in range(0, numofruns+1):
        top12t.append({})
        top12.append({})
    for x in range(1, numofruns+1):
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for k in range(mink, maxk+1):
            top12t[x][k] = list(sorted(kmercount[x][k].items(), key=lambda x: x[1], reverse=True))[:12]
            top12[x][k] = [j[0] for j in top12t[x][k]]

top12maker(numofruns)
#print(kmercount)
#print(top12)

top6all = []

def top6(run, k):
    top6us = []
    top6ts = []
    top6ds = {}
    top6s = []
    topvalpos = 0
    while len(top6us) <= 11:
        next = top12[run][k][topvalpos]
        nkmer = hash2kmer(next, k)
        nrkmer = revComp(nkmer)
        nhrkmer = kmer2hash(nrkmer)
        if next not in top6us:
            top6us.append(next)
            #print("k "+str(k)+" len "+str(len(top6s))+" kmer "+str(next))
        if next == nhrkmer:
            if len(top6us) <= 11:
                top6us.append(nhrkmer)
        try:
            if kmercount[run][k][nhrkmer] > 0:
                if nhrkmer not in top6us:
                    if len(top6us) <= 11:
                        top6us.append(nhrkmer)
                        #print("k "+str(k)+" len "+str(len(top6s))+" kmer "+str(nhrkmer))
        except:
            if len(top6us) <= 11:
                top6us.append(next)
        topvalpos += 1

    for i in range(len(top6us)):
        if i//2 == (i+1)//2:
            top6ts.append((top6us[i], top6us[i+1]))

    for x in top6ts:
        try:
            first = kmercount[run][k][x[0]]
        except:
            first = 0
        if x[0] != x[1]:
            try:
                second = kmercount[run][k][x[1]]
            except:
                second = 0
        if x[0] == x[1]:
            second = 0
        top6ds[x] = first+second

    sort = sorted(top6ds.items(), key=operator.itemgetter(1), reverse=True)

    for z in range(len(sort)):
        top6s.append(sort[z][0][0])
        top6s.append(sort[z][0][1])

    return top6s




def top6maker(numofruns):
    for _ in range(0, numofruns+1):
        top6all.append({})
    for x in range(1, numofruns+1):
        #mink = int(inputlist[((x-1)*5)+7])
        #maxk = int(inputlist[((x-1)*5)+8])
        for k in range(mink, maxk+1):
            top6all[x][k] = top6(x, k)

top6maker(numofruns)
#print("top6all")
#print(top6all)

totaldict = []

def findalltotals():
    global totaldict
    for _ in range(numofruns+1):
        totaldict.append({})
    for runnum in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            totaldict[runnum][k] = sum(kmercount[runnum][k].values())

findalltotals()

#print("Total")
#print(totaldict)

#print(kmercombinations)

#print("kmercount")
#print(kmercount)
#print("hamlist")
#print(hamlist)
#print("hamdict")
#print(hamdict)
#print("pwm")
#print(pwm)

#print(os.getcwd())

def organisefigures():
    cwd = str(os.getcwd())
    try:
        os.makedirs(cwd+"/figures/"+str(identifier))
        os.makedirs(cwd+"/figures/"+str(identifier)+"/logos")
        os.makedirs(cwd+"/figures/"+str(identifier)+"/hamming_distance")
        os.makedirs(cwd+"/figures/"+str(identifier)+"/position_bias")
        os.makedirs(cwd+"/figures/"+str(identifier)+"/kmer_frequency")
    except:
        pass

organisefigures()


print('Generating Logos')
import WebLogoMod

print('Generating Hamming distance figures')
import HamDistFig

print('Generating Position bias figures')
import PositionBias

if numofruns > 1:
    print('Generating Kmer frequency figures')
    import KmerFreq
