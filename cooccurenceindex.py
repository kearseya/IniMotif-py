from Bio import SeqIO
import numpy as np

file = str(input("File name with path: "))

sequence1 = str(input("Frist motif: "))
sequence2 = str(input("Second motif: "))
filterchoice = str(input("Noise filter? "))
if filterchoice in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
    filter = True
    gaplength = int(input("Gap size: "))
else:
    filter = False

seq1len = len(sequence1)
seq2len = len(sequence2)


forwardandreverse = 0
forwardonly = 0
reverseonly = 0
gaplist = []


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

def cooccurence(FileName):
    fileopen = open(FileName, "r")
    firstline = fileopen.readline()
    firstline = firstline.strip()
    if firstline.startswith(">"):
        filetype = "fasta"
    if firstline.startswith("@"):
        filetype = "fastq"
    seq1 = kmer2hash(sequence1)
    seq2 = kmer2hash(sequence2)
    global forwardandreverse
    global forwardonly
    global reverseonly
    if seq1len == seq2len:
        k = seq1len
        for sequence in SeqIO.parse(FileName, filetype):
            line = str(sequence.seq)
            done = set()
            noisefilter = {"F": [],"R": []}
            for x in range(0,((len(line)+1)-k)):
                kmers = str(line[x:x+k])
                try:
                    hkmer = kmer2hash(kmers)
                    done.add(hkmer)
                    if hkmer == seq1:
                        noisefilter["F"].append(x)
                    if hkmers == seq2:
                        noisefilter["R"].append(x)
                except:
                    continue
                if x == (len(line)-k):
                    if filter = True:
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
                            if k < mindiff <= (k+gaplength):
                                forwardandreverse += 1
                                gaplist.append(mindiff-k)
                        elif seq1 in done and seq2 not in done:
                            forwardonly += 1
                        elif seq2 in done and seq1 not in done:
                            reverseonly += 1
                    else:
                        if seq1 & seq2 in done:
                            forwardandreverse += 1
                        elif seq1 in done and seq2 not in done:
                            forwardonly += 1
                        elif seq2 in done and seq1 not in done:
                            reverseonly += 1


cooccurence(file)


def cooccurenceindex():
    cooccurence = (forwardandreverse/(forwardandreverse+forwardonly+reverseonly))
    return cooccurence

cooccurenceindex = cooccurenceindex()

print("Cooccurence index:")
print(cooccurenceindex)

def meangapcalc():
    if filter = False:
        meangap = "Not calculated"
    if filter = True:
        if len(gaplist) > 0:
            meangap = int(round((sum(gaplist)/len(gaplist)),0))
        else:
            meangap = 0
    return meangap

meangap = meangapcalc()

def stdevcalc():
    if filter = False:
        stdev = "Not calculated"
    if filter = True:
        if len(gaplist) > 0:
            stdev = statistics.stdev(gaplist)
        else:
            stdev = 0
    return stdev

stdev = stdevcalc()

print("Mean gap:")
print(meangap)
print("Standard deviation:")
print(stdev)
