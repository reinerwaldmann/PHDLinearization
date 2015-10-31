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


# Надлежит получить все данные о распределении каждой ячейки матрицы. SO:
# самое простое - получить набор гистограмм для каждого столбца, на столбец по папке.
# далее, критерий нормальности распределения для половины распределения
# (лол, самое простое - к тестируемому распределению добавить его вторую половинку, взятую с отрицательным знаком)

def makeListOfListPerIteration(listoflistoflistofdiffs):
    """
    takes a list of cases
    each case is a list of lists of end differences beginning from the first

    we make a list of list, like a end-difference matrix, where
     each item is a list.
    """
    #print ("number of types of differences", len(listoflistoflistofdiffs[0]))
    res= []
    for difflistnum in range (len(listoflistoflistofdiffs[0])): # number of enddiffs, actually fixed and is 4
        sq=[]
        for diffnum in range (len(listoflistoflistofdiffs[0][difflistnum])): #  diffnum - number of enddiff
            #print ("number of difference",len(listoflistoflistofdiffs[0][difflistnum]))
            listOfDiffsSameIterationPerCases = []
            for case in range (len(listoflistoflistofdiffs)):
                try:

                    listOfDiffsSameIterationPerCases.append(listoflistoflistofdiffs[case][difflistnum][diffnum])
                except:
                    #print ('shit happens', case, difflistnum,diffnum)
                    pass


            sq.append(listOfDiffsSameIterationPerCases)
        res.append(sq)


 #   print (res)
    return res

def printListOfListPerIteration(olololist):
    """
    just prints it
    """
    # cycle per diff
    j=1
    for diff in olololist:
        print ('Number of difference level: {0}'.format(j))
        i=0
        for iteration in diff:
            print ('Difference {0} - {1}'.format(i, iteration))
            i+=1
        print (" ")
        j+=1

    description = "See table of end differences, actually: number of difference level: \n" \
                "is a number of difference level (from 1 to 4, may be more, but probably not needed) \n" \
                "then, number of difference: if there'd be difference number 0, then it'l be number of \n" \
                "iteration, now it's the number of rows occupied by difference. \n" \
                "The list appearing next to each difference - each value per case. \n" \
                "As one can see, some cases finished in less number of iterations, so some lists are shorter \n" \
                "then others."



    print (description)


def makeAllEndDiffs (datafile):
    with open(datafile, 'rb') as f:
        rrlist = pickle.load(f)
        rrlist_filtered = [case for case in rrlist if case.ll()<30]
        return [transform_iia_to_matrix(lst, 0) for lst in rrlist_filtered ]



def makeAverageListOutOfEndDifferencesList(listdiffs):
    """
    """
    import copy
    import statistics
    lst = copy.copy(listdiffs)


    for j in lst:
        for m in j:
            for t in m:
                t = statistics.mean(t)
    return lst

def makeDispOutOfEndDifferences(listdiffs):
    """
    """
    pass

def makeHistsOutOfEndDifferences(listdiffs):
    """

    """

    pass



def test():

    datafile = "resfiles/resdump.dat"
    a = makeAllEndDiffs(datafile)

#    [print (r) for r in a]

    b = makeListOfListPerIteration(a)
    printListOfListPerIteration(b)




if __name__ == '__main__':
    test()




