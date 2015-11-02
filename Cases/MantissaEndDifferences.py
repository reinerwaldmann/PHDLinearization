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
        #mat[i*2][0] = abs(iia[i]['b'][index])
        mat[i*2][0] = iia[i]['b'][index]



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
            #mat[row][col]=abs(mat[row+1][col-1]-mat[row-1][col-1])
            mat[row][col]=mat[row+1][col-1]-mat[row-1][col-1]




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

def printListOfListPerIteration(olololist, file=None):
    """
    just prints it
    """
    # cycle per diff
    j=1
    for diff in olololist:
        print ('Number of difference level: {0}'.format(j), file=file)
        i=0
        for iteration in diff:
            print ('Difference {0} - {1}'.format(i, iteration), file=file)
            i+=1
        print (" ", file)
        j+=1

    description = "See table of end differences, actually: number of difference level: \n" \
                "is a number of difference level (from 1 to 4, may be more, but probably not needed) \n" \
                "then, number of difference: if there'd be difference number 0, then it'l be number of \n" \
                "iteration, now it's the number of rows occupied by difference. \n" \
                "The list appearing next to each difference - each value per case. \n" \
                "As one can see, some cases finished in less number of iterations, so some lists are shorter \n" \
                "then others."



    print (description, file=file)


def makeAllEndDiffs (datafile, index):
    with open(datafile, 'rb') as f:
        rrlist = pickle.load(f)
        rrlist_filtered = [case for case in rrlist if case.ll()<30]
        return [transform_iia_to_matrix(lst, index) for lst in rrlist_filtered ]

STAT_MEAN=0
STAT_VAR=1
STAT_SIGMA=2
STAT_SKEW=3
STAT_KURTOSIS=4



def makeCharactOutOfEndDifferencesList(listdiffs, type=STAT_MEAN):
    """
    type=0, mean
    1, variance
    2, sigma
    """
    import copy
    import statistics
    import scipy.stats as scps


    lst = copy.deepcopy(listdiffs)

    func={0:statistics.mean, 1:statistics.variance, 2:lambda x: math.sqrt(statistics.variance(x)), 3:scps.skewtest, 4:scps.kurtosistest  }

    for j in range (len(lst)):
        for k in range (len(lst[j])):
            lst[j][k] = func[type](lst[j][k])

    return lst





def makeHistsOutOfEndDifferences(listdiffs, mainfolder, limy=200):
    """
    limy - end value oy
    """
    import os, os.path

    import matplotlib.pyplot as plt
    # os.makedirs(path[, mode])¶

    for j in range (len(listdiffs)):  #для каждого уровня разности
        fldpath = os.path.join(os.path.dirname(mainfolder), "diflevel_{0}".format(j)) #  получили путь к папке
        os.makedirs(fldpath)

        for m in range(len(listdiffs[j])):
            fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
            # x1,x2,y1,y2 = plt.axis()
            # plt.axis((x1,x2,y1,limy))

            #   axes = plt.gca()
            #axes.set_ylim([0,200])


            #ax.hist( [x for x in listdiffs[j][m] if x<.003] ,30)
            ax.hist( [x for x in listdiffs[j][m] ] ,30)



            #это у нас графики значений компонента вектора по итерациям.
            imgpath = os.path.join(fldpath, "img_{0}.png".format(m))
            fig.savefig(imgpath)
            plt.close(fig)






from PyQt5.QtWidgets import QApplication, QWidget,  QMainWindow
from PyQt5 import uic


from uifiles.Mantissa2UI.m2mainwindow import Ui_MainWindow




class MainWindow(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.push1.clicked.connect(self.handle_button_1)


    def handle_button_1(self):
        self.textConsole.append("Button1 pressed")



def test():

    import os.path

    datafile = "resfiles/resdump205_DISP.dat"
    #datafile = "resfiles/resdump205.dat"

    #folder = "resfiles/hist_all_end_diffs_limhists/"

    folder = "resfiles/hists_all_end_diffs/"


    b = makeAllEndDiffs(datafile,2) # 0 - index of a variable, this case b0


    a = makeListOfListPerIteration(b)




    #
    # with open(os.path.join(folder,'average.txt'), 'wt') as f:
    #     printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_MEAN), f)
    #
    # with open(os.path.join(folder,'variance.txt'), 'wt') as f:
    #         printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_VAR), f)
    #
    # with open(os.path.join(folder,'sigma.txt'), 'wt') as f:
    #         printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_SIGMA), f)

    #makeHistsOutOfEndDifferences(a, folder, limy=200)
    #makeHistsOutOfEndDifferences(a, folder)

# NOW WE TEST FOR NORMALITY

    import scipy.stats as scps


    import random
    random.seed()

    normlist = [random.normalvariate(10, 2) for x in range(200)]

    print (scps.skewtest(normlist))
    print (scps.kurtosistest(normlist))
    print (scps.normaltest(normlist))




    exit(0)



    with open(os.path.join(folder,'skew.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_SKEW), f) #  critical 0.34  http://mvpprograms.com/help/mvpstats/distributions/SkewnessCriticalValues

    with open(os.path.join(folder,'kurtosis.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_KURTOSIS), f) #critical low -0.54	high 0.79  http://mvpprograms.com/help/mvpstats/distributions/KurtosisCriticalValues


# http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3591587/

    #https://onlinecourses.science.psu.edu/statprogram/node/138







# we must test got distribution for normality
# http://mvpprograms.com/help/mvpstats/distributions/KurtosisCriticalValues

# http://mvpprograms.com/help/mvpstats/contents

#http://mvpprograms.com/help/mvpstats/distributions/SkewnessCriticalValues








    # app = QApplication(sys.argv)
    # w = MainWindow()
    # w.show()
    # sys.exit(app.exec_())





#    [print (r) for r in a]

    #printListOfListPerIteration(makeAverageListOutOfEndDifferencesList(b))







if __name__ == '__main__':
    test()




