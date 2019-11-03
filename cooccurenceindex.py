from Bio import SeqIO
import numpy as np

file = str(input("File name with path: "))

sequence1 = str(input("Frist motif: "))
sequence2 = str(input("Second motif: "))

seq1len = len(sequence1)
seq2len = len(sequence2)


forwardandreverse = 0
forwardonly = 0
reverseonly = 0


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
            for x in range(0,((len(line)+1)-k)):
                kmers = str(line[x:x+k])
                try:
                    done.add(kmer2hash(kmers))
                except:
                    continue
                if x == (len(line)-k):
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

print(cooccurenceindex)
