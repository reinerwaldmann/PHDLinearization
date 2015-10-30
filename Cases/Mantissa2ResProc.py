__author__ = 'vasilev_is'

import pickle

import matplotlib.pyplot as plt
from Cases.CasesUtilStuff import IterationInfoAcceptor
import numpy as np
import math

import os, sys

def norm1 (b1,b2):
    absdif = [math.fabs(b11-b22) for b11,b22 in zip(b1,b2)]
    return max(absdif)

def norm2 (b1,b2):
    return np.linalg.norm(np.array(b1)-np.array(b2))

def norm3 (b1,b2):
    absdif = [math.fabs(b11-b22) for b11,b22 in zip(b1,b2)]
    return sum(absdif)



datafile = 'resfiles/resdump205_DISP.dat'

with open(datafile, 'rb') as f:
    rrlist = pickle.load(f)

#фильтруем случаи, где за 30 итераций не сошлось
rrlist_filtered = [case for case in rrlist if case.ll()<30]

print ('total', len(rrlist))

print ('dif', len(rrlist)-len(rrlist_filtered))


folder = "resfiles/chists_spread/chists_spread_n/"



for coeff in range (0,3):

    os.mkdir(folder+"changes_hists_b{0}".format(coeff))

    #для каждой итерации создаём гистограмму рассеяния изменений
    for it in range (10):
        #it=1
        dflist=[]  # список изменений на первой итерации по кейсам
        for case in rrlist_filtered:
            #dflist.append(norm3(case[0]['b'], case[1]['b']))
            try:
                dflist.append(case[it+1]['b'][coeff]- case[it]['b'][coeff])
            except:
                pass

        fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
        ax.hist(dflist,30)

        #это у нас графики значений компонента вектора по итерациям.
        fig.savefig (folder+'changes_hists_b{1}/img_{0}.png'.format(it, coeff))
        plt.close(fig)




