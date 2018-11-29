import seaborn
import matplotlib.pyplot as plt
plt.style.use('seaborn-ticks')
from matplotlib import transforms
import matplotlib.patheffects
import numpy as np

from test import pwm
from test import numofruns
from test import mink
from test import maxk



logoform = []


def logopos(a,t,c,g):
    pos =  [('T', t), ('G', g), ('A', a), ('C', c)]
    pos.sort(key = lambda x: x[1], reverse = False)
    return pos



def kmerpwm(runnum, k):
    kmerpwm = []
    for i in range(1,k+1):
        kmerpwm.append(logopos(pwm[runnum][k][i]["T"],pwm[runnum][k][i]["G"],pwm[runnum][k][i]["A"],pwm[runnum][k][i]["C"]))
    return kmerpwm



def allmaker(numofruns, mink, maxk):
    for z in range(1, numofruns+1):
        for j in range(mink,maxk+1):
            logoform.append(kmerpwm(z,j))

allmaker(numofruns, mink, maxk)

print(logoform)




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

def draw_logo(all_scores):
    fig = plt.figure()
    fig.set_size_inches(len(all_scores),2.5)
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


    seaborn.despine(ax=ax, offset=30, trim=True)
    ax.set_xticklabels(range(1,len(all_scores)+1), rotation=90)
    ax.set_yticklabels(np.arange(0,3,1))
    plt.show()




def logoprinter():
    for i in logoform:
        draw_logo(i)


logoprinter()
