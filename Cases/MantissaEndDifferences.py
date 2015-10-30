__author__ = 'vasilev_is'

import pickle

import matplotlib.pyplot as plt
from Cases.CasesUtilStuff import IterationInfoAcceptor
import numpy as np
import math

import os, sys

np.set_printoptions(precision=5, threshold=None, edgeitems=None, linewidth=300, suppress=None, nanstr=None, infstr=None, formatter=None)

def tt():

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


def transform_iia_to_matrix (iia, index):
    """

    :param iia:
    :param index: b is a vector, so changes a vector too. 'Cause of this we choose only one value of this vector, specified by index
    :return: np matrix, where 1st column filled with values on each iteration,
    not on each row, but separated with zero rows
    """

    mat = np.zeros((iia.ll()*2,5)) #numrows, numcolumns


    for i in range (iia.ll()):
        mat[i*2][ 0] = abs(iia[i]['b'][0])

    return makeEndDifferencesOnMatrix(mat,1)

def makeEndDifferencesOnMatrix(mat, squeeze=False):
    """
    fills matrix with end differences
    :param mat:
    :return:
    """
#     if x & 1:
#        return 'odd'
    # else:
    #    return 'even'

    # ROW-COLUMN!!!!!!!!

    for col in range (1,mat.shape[1]):
        #start=1 if col&1 else 0;
        start = col
        for row in range (start, mat.shape[0]-1-col, 2):
            mat[row][col]=abs(mat[row+1][col-1]-mat[row-1][col-1])

    if (squeeze):
        """
        Вернуть список списков конечных разностей
        """
        return [ [mat[row][col] for row in range (col, mat.shape[0]-1-col,2) ]  for col in range(1, mat.shape[1]) ]
    return mat



def test():

    datafile = "resfiles/resdump205_DISP.dat"

    with open(datafile, 'rb') as f:
        rrlist = pickle.load(f)
    rrlist_filtered = [case for case in rrlist if case.ll()<30]



    print (transform_iia_to_matrix(rrlist_filtered[0], 0))

if __name__ == '__main__':
    test()




