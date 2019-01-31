#%matplotlib inline
from KmerKounter import identifier
import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
import numpy as np
import math
#import operator

from KmerKounter import hamlist
from KmerKounter import pwm
from KmerKounter import numofruns
from KmerKounter import mink
from KmerKounter import maxk


def En(n):
    Enn = ((1/math.log(2))*((4-1)/(2*n)))
    return Enn

def Shannon(runnum, k, pos):
    bases = ["A", "T", "C", "G"]
    Hi = 0
    for b in bases:
        f = pwm[runnum][k][pos][b]
        Hi += -(f*math.log2(f))
    return Hi

def Ri(runnum, k, i):
    n = len(hamlist[runnum][k])
    Rii = math.copysign((math.log2(4)-(Shannon(runnum, k, i) + En(n))), 1)
    return Rii


logoform = []

def dictinit():
    for r in range(0, numofruns+1):
        logoform.append({})
        for k in range(mink, maxk+1):
            logoform[r].update({k:[]})
dictinit()


def logopos(a,t,c,g):
    upos =  [('A', a), ('T', t), ('C', c), ('G', g)]
    pos = sorted(upos, key=lambda x:x[1])
    return pos


def kmerpwm(runnum, k):
    kmerpwm = []
    for i in range(1,k+1):
        kmerpwm.append(logopos(pwm[runnum][k][i]["A"]*Ri(runnum, k, i),pwm[runnum][k][i]["T"]*Ri(runnum, k, i),pwm[runnum][k][i]["C"]*Ri(runnum, k, i),pwm[runnum][k][i]["G"]*Ri(runnum, k, i)))
    return kmerpwm





def allmaker(numofruns, mink, maxk):
    for z in range(1, numofruns+1):
        for j in range(mink,maxk+1):
            logoform[z][j].append(kmerpwm(z,j))

allmaker(numofruns, mink, maxk)




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

def draw_logo(all_scores, run, k):
    fig = plt.figure()
    fig.set_size_inches(len(all_scores)+1,2.5)
    ax = fig.add_subplot(111)
    ax.set_xticks(range(len(all_scores)))

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
    ax.sey_ylabel("bits")


    seaborn.despine(ax=ax, offset=30, trim=True)
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,3,1))
    #plt.show()
    plt.savefig('figures/logo_'+str(identifier)+"_"+str(run)+"_"+str(k), dpi=600)
    plt.close()




def logoprinter():
    for r in range(1, numofruns+1):
        for i in logoform[r]:
            for j in range(0, len(logoform[r][i])):
                draw_logo(logoform[r][i][j], r, i)


logoprinter()
