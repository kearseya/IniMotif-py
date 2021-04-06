import numpy as np
import random

bases = ["A", "T", "C", "G"]
motif = "ATCGTATACC"

def stringmaker():
    string = ""
    for _ in range(100):
        string += bases[round(random.uniform(0,3))]
    return(string)

def motifinset(string):
    position = round(random.uniform(0,90))
    insert = string.replace(string[position:position+10], motif)
    return(insert)


def fastamaker():
    with open("simdata", "w") as f:
        for _ in range(100000):
            header = ">header"
            dna = motifinset(stringmaker())
            f.write(header+'\n')
            f.write(dna+'\n')

fastamaker()
