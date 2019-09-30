#%matplotlib inline
from KmerKounter import identifier
import seaborn
import matplotlib
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
import numpy as np
import math
#import operator

from KmerKounter import revnuc
from KmerKounter import revComp
#from KmerKounter import hamlist
#from KmerKounter import pwm
from KmerKounter import numofruns
from KmerKounter import mink
from KmerKounter import maxk
from KmerKounter import startround

from KmerKounter import top6all
from KmerKounter import kmercount
from KmerKounter import totaldict

from KmerKounter import inputlist

matplotlib.use('AGG')

"""
from KmerKounter import mink
from KmerKounter import maxk
"""

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



def inputtype():
    global nmotifs
    global type
    global allowham
    try:
        from GUI import nmotifs
        from GUI import logotype
        #if len(logotype) > 1:
        type = logotype[0]
        from GUI import allowham
    except:
        nmotifs = int(input("Number of motifs: "))
        type = input("Bits (b) or Frequency (f)?: ")
        allowham = int(input("Number of motifs: "))

inputtype()

#print("PWM")
#print(pwm)

hamlist = []
rhamlist = []
removelist = []
hamdict = []
rhamdict = []
countdict = []
pwm = []
logoform = []

def dictinit(numofruns):
    for _ in range(numofruns+1):
        #kmercount.append({})
        hamlist.append({})
        rhamlist.append({})
        removelist.append({})
        hamdict.append({})
        rhamdict.append({})
        countdict.append({})
        pwm.append({})
        logoform.append({})

dictinit(numofruns)


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def removelistinit():
    for runnum in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            removelist[runnum][k] = {}
            for n in range(1, nmotifs+1):
                removelist[runnum][k][n] = set()

removelistinit()


def listhammer():
    for runnum in range(1, numofruns+1):
        for i in kmercount[runnum]:
            hamlist[runnum][i] = {}
            rhamlist[runnum][i] = {}
            #removelist[runnum][i] = {}
            for n in range(1, nmotifs+1):
                hamlist[runnum][i][n] = []
                rhamlist[runnum][i][n] = []
                #removelist[runnum][i][n] = set()
                #hconsensus = (list(kmercount[runnum][i].keys())[0])
                if n == 1:
                    hconsensus = max(kmercount[runnum][i], key=lambda key: kmercount[runnum][i][key])
                else:
                    temptop6 = top6all[runnum][i].copy()
                    for x in temptop6:
                        if x in removelist[runnum][i][n]:
                            temptop6.remove(x)
                            temptop6.remove(kmer2hash(revComp(hash2kmer(x, i))))
                    hconsensus = max(temptop6, key=lambda key: kmercount[runnum][i][key])
                consensus = hash2kmer(hconsensus, i)
                for x in list(kmercount[runnum][i].keys()):
                    values = hash2kmer(x,i)
                    if hamming_distance(consensus, values) <= allowham:
                        if x not in removelist[runnum][i][n]:
                            rvalue = revComp(values)
                            rx = kmer2hash(rvalue)
                            try:
                                if kmercount[runnum][i][x] >= kmercount[runnum][i][rx]:
                                    hamlist[runnum][i][n].append(x)
                                    for j in range(n+1, nmotifs+1):
                                        removelist[runnum][i][j].add(x)
                                        removelist[runnum][i][j].add(rx)
                                    if x != rx:
                                        rhamlist[runnum][i][n].append(rx)
                                        if n != nmotifs:
                                            for j in range(n+1, nmotifs+1):
                                                removelist[runnum][i][j].add(x)
                                                removelist[runnum][i][j].add(rx)
                            except:
                                hamlist[runnum][i][n].append(x)
                                for j in range(n+1, nmotifs+1):
                                    removelist[runnum][i][j].add(x)
                                    removelist[runnum][i][j].add(rx)

listhammer()

"""
def visualise():
    visibleremove = []
    for _ in range(0, numofruns+1):
        visibleremove.append({})
    for r in range(1, numofruns+1):
        visibleremove[r] = {}
        for k in range(mink, maxk+1):
            visibleremove[r][k] = {}
            for n in range(1, nmotifs+1):
                visibleremove[r][k][n] = []
                for x in removelist[r][k][n]:
                    visibleremove[r][k][n].append(hash2kmer(x, k))
    return visibleremove
"""
#print(visualise())

#print(removelist)

"""
else:
    po = 0
    for p, x in enumerate(top6all[runnum][i], 1):
        if p > po:
            if x not in removelist[runnum][i][n]:
                hconsensus = x
                po = p
"""

#print(hamlist)
#print(rhamlist)
#print("remove")
#print(removelist)

def dicthammer():
    global countdict
    for runnum in range(1, numofruns+1):
        for k in range(mink, maxk+1):
            hamdict[runnum][k] = {}
            rhamdict[runnum][k] = {}
            countdict[runnum][k] = {}
            for n in range(1, nmotifs+1):
                hamdict[runnum][k][n] = {}
                rhamdict[runnum][k][n] = {}
                countdict[runnum][k][n] = 0
                frac = { z : kmercount[runnum][k][z] for z in hamlist[runnum][k][n] }
                rfrac = { z : kmercount[runnum][k][z] for z in rhamlist[runnum][k][n] }
                combine = {**frac, **rfrac}
                #M = sum(kmercount[runnum][k].values())
                M = sum(combine.values())
                frac = { z : (kmercount[runnum][k][z]/M) for z in frac }
                rfrac = { z : (kmercount[runnum][k][z]/M) for z in rfrac }
                hamdict[runnum][k][n].update(frac)
                rhamdict[runnum][k][n].update(rfrac)
                for x in hamlist[runnum][k][n]:
                    countdict[runnum][k][n] += kmercount[runnum][k][x]
                for y in rhamlist[runnum][k][n]:
                    countdict[runnum][k][n] += kmercount[runnum][k][y]

dicthammer()

#print("Count")
#print(countdict)

def startpwm():
    for runnum in range(1, numofruns+1):
        for i in hamdict[runnum]:
            pwm[runnum][i] = {}
            logoform[runnum][i] = {}
            for n in range(1, nmotifs+1):
                pwm[runnum][i][n] = {}
                logoform[runnum][i][n] = []
                for j in range(1,i+1):
                    pwm[runnum][i][n][j] = {"A":0}
                    pwm[runnum][i][n][j].update({"C":0})
                    pwm[runnum][i][n][j].update({"G":0})
                    pwm[runnum][i][n][j].update({"T":0})

startpwm()

def pwmmaker():
    for runnum in range(1, numofruns+1):
        for i in hamdict[runnum]:
            for n in range(1, nmotifs+1):
                for x in hamdict[runnum][i][n]:
                    kmer = hash2kmer((x),i)
                    for j in range(1,i+1):
                        if kmer[j-1] == "A":
                            pwm[runnum][i][n][j]["A"] += hamdict[runnum][i][n][x]
                        elif kmer[j-1] == "C":
                            pwm[runnum][i][n][j]["C"] += hamdict[runnum][i][n][x]
                        elif kmer[j-1] == "T":
                            pwm[runnum][i][n][j]["T"] += hamdict[runnum][i][n][x]
                        elif kmer[j-1] == "G":
                            pwm[runnum][i][n][j]["G"] += hamdict[runnum][i][n][x]
                for x in rhamdict[runnum][i][n]:
                    kmer = revComp(hash2kmer((x),i))
                    for j in range(1,i+1):
                        if kmer[j-1] == "A":
                            pwm[runnum][i][n][j]["A"] += rhamdict[runnum][i][n][x]
                        elif kmer[j-1] == "C":
                            pwm[runnum][i][n][j]["C"] += rhamdict[runnum][i][n][x]
                        elif kmer[j-1] == "T":
                            pwm[runnum][i][n][j]["T"] += rhamdict[runnum][i][n][x]
                        elif kmer[j-1] == "G":
                            pwm[runnum][i][n][j]["G"] += rhamdict[runnum][i][n][x]

pwmmaker()

#print(pwm)

def En(num):
    Enn = ((1/math.log(2))*((4-1)/(2*num)))
    return Enn

def Shannon(runnum, k, pos, n):
    bases = ["A", "T", "C", "G"]
    Hi = 0
    for b in bases:
        f = pwm[runnum][k][n][pos][b]
        try:
            Hi += -(f*math.log2(f))
        except:
            continue
    return Hi

def Ri(runnum, k, i, n):
    num = len(hamlist[runnum][k][n])+len(rhamlist[runnum][k][n])
    Rii = math.copysign((math.log2(4)-(Shannon(runnum, k, i, n) + En(num))), 1)
    return Rii



"""
def dictinit():
    for r in range(0, numofruns+1):
        logoform.append({})
        for k in range(mink, maxk+1):
            logoform[r][k] = {}
            #logoform[r].update({k:[]})
            for n in range(1, nmotifs+1):
                logoform[r][k][n] = []
dictinit()
"""

def logopos(a,t,c,g):
    upos =  [('A', a), ('T', t), ('C', c), ('G', g)]
    pos = sorted(upos, key=lambda x:x[1])
    return pos

#print("PWM")
#print(pwm)
#print("LOGO")
#print(logoform)


def kmerpwm(runnum, k, n):
    kmerpwm = []
    for i in range(1,k+1):
        if type == "f":
            kmerpwm.append(logopos(pwm[runnum][k][n][i]["A"],pwm[runnum][k][n][i]["T"],pwm[runnum][k][n][i]["C"],pwm[runnum][k][n][i]["G"]))
        if type == "b":
            kmerpwm.append(logopos(pwm[runnum][k][n][i]["A"]*Ri(runnum, k, i, n),pwm[runnum][k][n][i]["T"]*Ri(runnum, k, i, n),pwm[runnum][k][n][i]["C"]*Ri(runnum, k, i, n),pwm[runnum][k][n][i]["G"]*Ri(runnum, k, i, n)))
    return kmerpwm



def allmaker(numofruns):
    for z in range(1, numofruns+1):
        #mink = int(inputlist[((z-1)*5)+7])
        #maxk = int(inputlist[((z-1)*5)+8])
        for j in range(mink,maxk+1):
            for n in range(1, (nmotifs+1)):
                logoform[z][j][n].append(kmerpwm(z,j,n))

allmaker(numofruns)


#print("Logoform")
#print(logoform)


COLOR_SCHEME = {'G': 'orange',
                'A': 'red',
                'C': 'blue',
                'T': 'darkgreen'}
BASES = list(COLOR_SCHEME.keys())


class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def draw_logo(all_scores, run, k, n):
    fig = plt.figure()
    fig.set_size_inches(k, 2)
    ax = fig.add_subplot(111)
    #ax.set_xticks(range(k))

    xshift = 0
    trans_offset = transforms.offset_copy(ax.transAxes,
                                      fig=fig,
                                      x=0,
                                      y=0,
                                      units='points')


    for scores in all_scores:
        yshift = 0
        for base, score in scores:
            txt = ax.text(0,
                          0,
                          base,
                          transform = trans_offset,
                          fontsize = 80,
                          color = COLOR_SCHEME[base],
                          weight = 'bold',
                          ha = 'center',
                          family = 'sans-serif'
                          )
            txt.set_clip_on(False)
            txt.set_path_effects([Scale(1.0, score)])
            fig.canvas.draw()
            window_ext = txt.get_window_extent(txt._renderer)
            yshift = window_ext.height*score
            trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='points')
        xshift += window_ext.width
        trans_offset = transforms.offset_copy(ax.transAxes, fig=fig, x=xshift, units='points')


    ax.set_yticks(range(0,3))
    if type == "f":
        ax.set_ylabel("frequency")
    if type == "b":
        ax.set_ylabel("bits")
    #ax.set_title("Uses "+str(round((countdict[run][k][n]/totaldict[run][k])*100, 2))+"%")
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    plt.tick_params(bottom = False, labelbottom = False)
    seaborn.despine(ax=ax, offset=30, trim=False, bottom=True)
    #ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,3,1))
    #plt.show()
    plt.savefig("figures/"+str(identifier)+"/logos/logo_"+str(identifier)+"_"+str(run+(startround-1))+"_"+str(k)+"_"+str(n), dpi=600, bbox_inches='tight')
    plt.close()



def logoprinter():
    for r in range(1, numofruns+1):
        for i in logoform[r]:
            for n in range(1, nmotifs+1):
                for j in range(0, len(logoform[r][i][n])):
                    draw_logo(logoform[r][i][n][j], r, i, n)


logoprinter()
