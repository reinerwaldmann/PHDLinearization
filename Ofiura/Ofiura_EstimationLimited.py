__author__ = 'vasilev_is'

import math, copy, os, platform


import numpy as np

import Ofiura.Ofiura_Qualitat as o_q
import Ofiura.Ofiura_general as o_g


"""
ПОстроена на использовании библиотеки numpy
"""




def makeSkInit (funcf, measdata, binit, c):
    Sk=None #или mpm.mpf(0.0)
    for point in measdata:
        fxbc=funcf(point['x'],binit,c)
        if fxbc is None:
            print ("ERROR: makeSkInit: funcf returned None")
            return None
        dif=point['y']-fxbc
        if Sk is None:
            Sk=np.dot(dif.T,dif)
        else:
            Sk+=np.dot(dif.T,dif)
    return Sk

def islinux ():
     if 'Windows' in platform.system():
         return False
     return True




def makeAinit (bstart, bend, Skbinit, binit):
    """
    расчёт ведётся по правилу "binit есть хорошее предположение в рамках области"
    :param bstart:
    :param bend:
    :param Skbinit:
    :param binit:
    :return:
    """


    A=np.zeros ( (len(binit), 2 ))  #,  dtype=np.float96 if islinux() else np.float64

    for i in range (len(binit)):
        #A[i][0]=A[i][1]=.001*(bend[i]-bstart[i])*Skbinit  #универсально
        A[i][0]=0.001*(binit[i]-bstart[i])*Skbinit
        A[i][1]=0.001*(bend[i]-binit[i])*Skbinit

    return A



def countSklims(A,b,bstart, bend):
    partone=parttwo=.0
    for i in range (len(b)):
        partone+=A[i][0]/(b[i]-bstart[i]) #UNSAFE: если вдруг, ну вдруг b=bstart или bend, то будет <strike> треш </strike> креш с делением на ноль.
        parttwo+=A[i][1]/(bend[i]-b[i])
    return partone+parttwo



def countN (A, b, bstart, bend):
    N=np.zeros (len(b), len(b))
    for j in range (len(b)):
        partone=parttwo=0
        for i in range (len(b)):
            partone+=2*A[i][0]/(b[i]-bstart[i])**3
            parttwo+=2*A[i][1]/(bend[i]-b[i])**3
        N[j][j]=parttwo+partone
    return N

