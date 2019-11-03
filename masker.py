import re
#from Bio import SeqIO
import os
import shutil
import fileinput

from itertools import chain, combinations, product

revnuc = {'A':'T','T':'A','G':'C','C':'G','N':'N'}

def revComp(seq):
    rev = ''
    for i in range(len(seq) - 1,-1,-1):
        rev += revnuc[seq[i]]
    return rev

def hamming_circle(s, n, alphabet):
    for positions in combinations(range(len(s)), n):
        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)
            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]
                else:
                    cousin[p] = alphabet[r]
            yield ''.join(cousin)

def hamming_ball(s, n, alphabet):
    allowedkmers = set(chain.from_iterable(hamming_circle(s, i, alphabet) for i in range(n + 1)))
    return allowedkmers

def myrepl(m):
    return ('N'*len(m.group(0)))

#pat = re.compile(f'({unit_str}){{{n_min_unit},{n_max_unit}}}')
#out_str = re.sub(pat, myrepl, instr)
max_rep = ''

file = str(input("File to be masked (with path): "))
outfilename =str(os.path.dirname(file))+"/masked_"+str(os.path.basename(file))

shutil.copy(file, str(os.path.dirname(file))+"/masked_"+str(os.path.basename(file)))

def masker():
    outfile = open(outfilename, "r")
    #outfile = open(str(os.path.dirname(file))+"/masked_"+str(os.path.basename(file)), "w")
    numberofmasks = int(input("Number of masks: "))
    for n in range(1, numberofmasks+1):
        typeofmask = str(input("Mask type (repeat/motif): "))
        if typeofmask in ["repeat", "r", "rep"]:
            unit = str(input("Unit string: "))
            revwant = str(input("Mask revcomp of unit?: "))
            min_rep = int(input("Minimum repeats: "))
            funit = {unit}
            pat = re.compile(f'({unit}){{{min_rep},{max_rep}}}')
            frunit = set()
            if revwant in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
                for i in funit:
                    frunit.add(i)
                    frunit.add(revComp(i))
                for unit in frunit:
                    pat = re.compile(f'({unit}){{{min_rep},{max_rep}}}')
                    for line in fileinput.input([outfilename], inplace=True):
                        print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))
            else:
                pat = re.compile(f'({unit}){{{min_rep},{max_rep}}}')
                for line in fileinput.input([outfilename], inplace=True):
                    print(line.replace(line, str.strip(re.sub(pat, myrepl, l)))
                    
        if typeofmask in ["motif", "m", "mot"]:
            unit = str(input("Unit string: "))
            min_rep = 1
            revwant = str(input("Mask revcomp of unit?: "))
            allowham = int(input("Hamming distance: "))
            maskmotifs = hamming_ball(unit, allowham, "ATGC")
            rmaskmotifs = set()
            if revwant in ["y", "Y", "yes", "Yes", "t", "true", "True"]:
                for i in maskmotifs:
                    rmaskmotifs.add(i)
                    rmaskmotifs.add(revComp(i))
                for kmer in rmaskmotifs:
                    pat = re.compile(f'({kmer}){{{min_rep},{max_rep}}}')
                    for line in fileinput.input([outfilename], inplace=True):
                        print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))
            else:
                for kmer in maskmotifs:
                    pat = re.compile(f'({kmer}){{{min_rep},{max_rep}}}')
                    for line in fileinput.input([outfilename], inplace=True):
                        print(line.replace(line, str.strip(re.sub(pat, myrepl, line))))

masker()
