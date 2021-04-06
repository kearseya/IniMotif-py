import itertools
#import csv
import operator

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.pyplot import imshow


filename = str(input("filename:"))

DNAbases = ["A", "T", "C", "G"]
all = {}
genecounts = {}
genetable = {}
krange = []
heatarrays = {}
maxcounts = {}
mincounts = {}


mink = int(input("Min K: "))
maxk = int(input("Max K: "))

mingl = int(input("Min gene length: "))
maxgl = int(input("Max gene length: "))



def initkrange(mink, maxk):
    for k in range(mink, maxk+1):
        print(k)
        krange.append(k)
initkrange(mink, maxk)

def initall():
    for k in krange:
        all[k] = [''.join(c) for c in itertools.product(DNAbases, repeat=k)]
initall()


def initgenetypes():
    for k in krange:
        genetable[k] = {}
        for s in all[k]:
            for e in all[k]:
                genetable[k][s+" "+e] = 0
initgenetypes()


def sequencefinder(filename):
    file = open(filename, 'r')
    igenome = "".join(line.strip() for line in file if line[0] not in ">")
    return igenome


def geneparser():
    genome = sequencefinder(filename)
    genomesize = len(genome)
    genes = 0
    print(genome[0:21])
    print("Total genome length: "+str(genomesize))
    for k in krange:
        genes = 0
        for s in all[k]:
            for e in all[k]:
                for x in range(0, genomesize-k):
                    kmer = genome[x:x+k]
                    if kmer == s:
                        #print(kmer)
                        for j in range(1, round(maxgl/k)):
                            nkmer = genome[x+(j*k):x+((j+1)*k)]
                            if nkmer in all[k]:
                                if nkmer not in [s, e]:
                                    #print("N "+nkmer)
                                    continue
                                if nkmer == e:
                                    if j*k >= mingl:
                                        #print("E "+nkmer)
                                        genes += 1
                                        genetable[k][kmer+" "+nkmer] += 1
                                        #print(genetable[k][kmer+" "+nkmer])
                                        break
                                if nkmer == s:
                                    #print("Restart")
                                    break
                            if nkmer not in all[k]:
                                #print("Invalid")
                                break
        genecounts[k] = genes
geneparser()

print(genecounts)
print(genetable)


def minfinder():
    for k in krange:
        minv = min(genetable[k], key=lambda key: genetable[k][key])
        mincounts[k] = genetable[k][minv]
minfinder()

def maxfinder():
    for k in krange:
        maxv = max(genetable[k], key=lambda key: genetable[k][key])
        maxcounts[k] = genetable[k][maxv]
maxfinder()


def heatstruct():
    for k in krange:
        heatarrays[k] = []
        for s in all[k]:
            row = []
            for e in all[k]:
                row.append(genetable[k][s+" "+e])
            heatarrays[k].append(row)
heatstruct()

#print(heatarrays)



def heatplotter():
    for k in krange:
        heatarraysnp = np.array(heatarrays)
        cmap = cm.Reds

        fig, ax = plt.subplots()
        im = plt.imshow(heatarrays[k], cmap=cmap)

        ax.set_xticks(np.arange(len(all[k])))
        ax.set_yticks(np.arange(len(all[k])))

        ax.set_xticklabels(all[k])
        ax.set_yticklabels(all[k])

        ax.set_ylabel("Start codons")
        ax.set_xlabel("End codons")

        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
            rotation_mode="anchor")

        cbar = ax.figure.colorbar(im, ax=ax, cmap=cmap)
        cbar.ax.set_ylabel("Counts", rotation=-90, va="bottom")


        ax.set_title("Start stop codon heatplot")
        fig.tight_layout()
        imshow(heatarrays[k], cmap=cmap)
        plt.show()

heatplotter()
