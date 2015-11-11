# -*- coding: utf-8 -*-

__author__ = 'vasilev_is'

import pickle
import statistics
import math
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np



# устанавливаем стиль fantasy по умолчанию
mpl.rcParams['font.family'] = 'fantasy'

# устанавливаем в стиль fantasy шрифт, который отображает кириллицу, например, Arial
mpl.rcParams['font.fantasy'] = 'Arial'


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
            print ('Difference {0} {1} values - {2}'.format(i,len(iteration) if hasattr(iteration,'__len__') else 0, iteration), file=file)
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
        #rrlist_filtered = [case for case in rrlist if case.ll()<30]
        return [transform_iia_to_matrix(lst, index) for lst in rrlist]

STAT_MEDIAN=-1
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

    func={-1:  statistics.median, 0:statistics.mean, 1:statistics.variance, 2:lambda x: math.sqrt(statistics.variance(x)), 3:scps.skewtest, 4:scps.kurtosistest  }

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
            fig, ax = plt.subplots( nrows=2, ncols=1 )  # create figure & 1 axis


            #ax[0].set_ylim([0,100])
            #ax[1].set_ylim([0,100])

            ax[0].hist( listdiffs[j][m], 30)
            mean = statistics.mean(listdiffs[j][m])
            median = statistics.median(listdiffs[j][m])

            rnd=10

            ax[0].set_title('Raw, mean={0}, median={1}'.format(round(mean, rnd), round(median, rnd)))


            cleaned = clean_exponent(listdiffs[j][m])


            ax[1].hist(cleaned ,30)
            mean = statistics.mean(cleaned)
            median = statistics.median(cleaned)

            number_of_deleted = abs(len(cleaned) - len(listdiffs[j][m]))
            percent_of_deleted = round(100*number_of_deleted/len(listdiffs[j][m]), 3)


            ax[1].set_title('Cleaned 3, mean={0}, medn={1}, del={2}/{3}%'.
                            format(round(mean,rnd), round(median,rnd), number_of_deleted, percent_of_deleted))


            #это у нас графики значений компонента вектора по итерациям.
            imgpath = os.path.join(fldpath, "img_{0}.png".format(m))
            fig.savefig(imgpath)
            plt.close(fig)

def makeHistsOutOfEndDifferences_simple(listdiffs, mainfolder, limy=200):
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


            #ax[0].set_ylim([0,100])
            #ax[1].set_ylim([0,100])

            ax.hist( listdiffs[j][m], 30)
            mean = statistics.mean(listdiffs[j][m])
            median = statistics.median(listdiffs[j][m])

            rnd=10

            ax.set_title('Raw, mean={0}, median={1}'.format(round(mean, rnd), round(median, rnd)))

            #это у нас графики значений компонента вектора по итерациям.
            imgpath = os.path.join(fldpath, "img_{0}.png".format(m))
            fig.savefig(imgpath)
            plt.close(fig)




def get_normalized_exp(x):
        """
        Получает экспоненту у нормализованного представления числа
        2 этапа: если оно имеет порядок 1 и выше, то срабатывает первая часть,
        до первого return (она уменьшает число последовательно)

        Иначе к делу подключается вторая

        """
        if x==0:
            return 0
            # or return False

        i=0
        while (int(x)):
            x/=10
            i+=1

        if i:
             return i


        while not (int(x)):
            x*=10
            i+=1

        return -1*i

def compare_exp (a,b, threshold=0):
    """
    returns
    1 (True), if a==b
    else 0 (False)
    Looking only at the exponent
    """

    return abs(get_normalized_exp(a)-get_normalized_exp(b))<=threshold


def list_to_exponents(listdiffs):
    """
    all values  in our diflist now will be exponents
    :param inlist:
    :return:
    """
    import copy

    lst = copy.deepcopy(listdiffs)
    #func={0:clean_Smirnoff, 1:clean_percentile}

    for j in range (len(lst)):
        for k in range (len(lst[j])):
            lst[j][k] = [get_normalized_exp(x) for x in lst[j][k]]

    return lst



def clean_MAD(inlist):
        median =  statistics.median(inlist)
        MAD = statistics.median([ abs(x-median)for x in inlist ])
        threshold = 3 # see pic on wiki from references above, this means like 3*sigma.
        # I suppose we can limit it even harder
        return [x for x in inlist if 0.6745*(abs(x-median)/MAD) < threshold ]

def clean_exponent(inlist):
    """
    cleans data thresh 2
    :param inlist:
    :return:
    """
    median = statistics.median(inlist)
    return [x for x in inlist if compare_exp(x,median, 0)]

def cleanData(listdiffs):
    """
    returns cleaned data.
    There are two strategies of cleaning:
        + clean data using Smirnoff criteria http://www.arhiuch.ru/lab4.html
        + clean data using simple method: reject 5-10%. Works this way: arrange values in a row sorted
            by their difference with mean value, then reject 5-10$ of the furthest.

    Information on the web:
    http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data
    "The problem with using percentile is that the points identified as outliers is a function of your sample size.
    There are a huge number of ways to test for outliers, and you should give some thought to how you classify them.
    Ideally, you should use a-priori information (e.g. "anything above/below this value is unrealistic because...")
    However, a common, not-too-unreasonable outlier test is to remove points based on their "median absolute deviation"."

     References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

        https://en.wikipedia.org/wiki/Standard_score#/media/File:Normal_distribution_and_scales.gif

        http://blog.caseystella.com/pyspark-openpayments-analysis-part-4.html

        http://stackoverflow.com/questions/22354094/pythonic-way-of-detecting-outliers-in-one-dimensional-observation-data


    :param listdiffs: input list
    :return: cleaned data list
    """
    import copy


    lst = copy.deepcopy(listdiffs)


    #func={0:clean_Smirnoff, 1:clean_percentile}
    func = clean_exponent

    for j in range (len(lst)):
        for k in range (len(lst[j])):
            lst[j][k] = func(lst[j][k])
    return lst


from PyQt5.QtWidgets import QMainWindow
from uifiles.Mantissa2UI.m2mainwindow import Ui_MainWindow

class MainWindow(QMainWindow, Ui_MainWindow):


    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.push1.clicked.connect(self.handle_button_1)


    def handle_button_1(self):
        self.textConsole.append("Button1 pressed")


def cleanFolder (folder):
      # cleaning res directory
    for subj in os.listdir(folder):
        if os.path.isdir(os.path.join(folder,subj)):
            cleanFolder(os.path.join(folder,subj))
            os.removedirs(os.path.join(folder,subj))
        else:
            os.remove(os.path.join(folder,subj))




def write_stat_characteristics(a, folder):
    with open(os.path.join(folder, 'median.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_MEDIAN),
                                    f)  # critical 0.34  http://mvpprograms.com/help/mvpstats/distributions/SkewnessCriticalValues
    with open(os.path.join(folder, 'mean.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_MEAN),
                                    f)  # critical low -0.54	high 0.79  http://mvpprograms.com/help/mvpstats/distributions/KurtosisCriticalValues
    with open(os.path.join(folder, 'skew.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_SKEW),
                                    f)  # critical 0.34  http://mvpprograms.com/help/mvpstats/distributions/SkewnessCriticalValues
    with open(os.path.join(folder, 'kurtosis.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_KURTOSIS),
                                    f)  # critical low -0.54	high 0.79  http://mvpprograms.com/help/mvpstats/distributions/KurtosisCriticalValues


    # make a row list
    with open(os.path.join(folder, 'raw.txt'), 'wt') as f:
        printListOfListPerIteration(a,
                                    f)  # critical low -0.54	high 0.79  http://mvpprograms.com/help/mvpstats/distributions/KurtosisCriticalValues
    with open(os.path.join(folder, 'average.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_MEAN), f)
    with open(os.path.join(folder, 'variance.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_VAR), f)
    with open(os.path.join(folder, 'sigma.txt'), 'wt') as f:
        printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a, STAT_SIGMA), f)


def makeBoxPlot(listdiff, mainfolder, index=0):
    """
    idea is to make boxplot for every iteration,
     so the trend will be seen.
    :param listdiff:
    :return:
    """
    # Create a figure instance
    fig = plt.figure(1, figsize=(9, 6))

    # Create an axes instance
    ax = fig.add_subplot(111)

    ax.set_title(u"Диаграмма разброса экспонент изменений по итерациям")

    # Create the boxplot
    #bp = ax.boxplot(listdiff[0])

    ## add patch_artist=True option to ax.boxplot()
    ## to get fill color
    bp = ax.boxplot(listdiff[0], patch_artist=True)

    ## change outline color, fill color and linewidth of the boxes
    for box in bp['boxes']:
        # change outline color
        box.set( color='#7570b3', linewidth=2)
        # change fill color
        box.set( facecolor = '#82D5FF' )

    ## change color and linewidth of the whiskers
    for whisker in bp['whiskers']:
        whisker.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the caps
    for cap in bp['caps']:
        cap.set(color='#7570b3', linewidth=2)

    ## change color and linewidth of the medians
    # for median in bp['medians']:
    #     median.set(color='#b2df8a', linewidth=2)

    ## change the style of fliers and their fill
    for flier in bp['fliers']:
        flier.set(marker='o', color='#e7298a', alpha=0.5)



    # Save the figure
    fig.savefig(os.path.join(mainfolder, "boxplot_1stdiff.png"), bbox_inches='tight')



def lstPerPercentil(lst, perc):
    import itertools
    lst.sort()
    keys = []
    groups=[]
    for k, g in itertools.groupby(lst):
        keys.append(k)
        groups.append(g)
        #print (k,g)

    print ('Порог', 'Число элементов > порога','Закрываемый процент', sep=',')
    superdict={}
    for k in keys:
        ll=len([x for x in lst if x>=k])
        #print ([x for x in lst if x>=k])
        print (k, ll , int(100*ll/len(lst)),sep=','  )
        superdict[int(100*ll/len(lst))]=k


    #print ()
    return superdict



def test():

    datafile = "resfiles/resdump205.dat"
    #datafile = "resfiles/resdump205.dat"
    #folder = "resfiles/hist_all_end_diffs_limhists/"

    folder = "resfiles/exponents/"


    with open(datafile, 'rb') as f:
        rl=pickle.load(f)

    for r in rl:
        for m in r:
            print (m)
    exit(0)


    #folder = "resfiles/hists_all_end_diffs/"


    b = makeAllEndDiffs(datafile,2) # second parameter - index of a variable
    a = makeListOfListPerIteration(b)
    a = list_to_exponents(a)

    i=0
    ddd={}
    for d in a[0]:
        i+=1
        print ('Iteration',i)
        ddd[i]=lstPerPercentil(d,90)


    print ('optimal')
    optperc=90

    for k in range (1,i+1):
        srtk = sorted(list(ddd[k].keys()))
        for f in srtk:
            if (f>=optperc):
                print(ddd[k][f])
                break









    # if (os.path.exists(folder)):
    #     cleanFolder(folder)
    # else:
    #     os.mkdir(folder)
    #
    # makeBoxPlot(a, folder)

#    write_stat_characteristics(a, folder)

 #   makeHistsOutOfEndDifferences_simple(a, folder)

    return (0)


# NOW WE TEST FOR NORMALITY
#
#     import scipy.stats as scps
#
#
#     import random
#     random.seed()
#
#     normlist = [random.normalvariate(10, 2) for x in range(200)]
#
#     print (scps.skewtest(normlist))
#     print (scps.kurtosistest(normlist))
#     print (scps.normaltest(normlist))
#     #kurtosis=эксцесс
#     #skew = ассиметричность
#
#
#
#
#
#     with open(os.path.join(folder,'skew.txt'), 'wt') as f:
#         printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_SKEW), f) #  critical 0.34  http://mvpprograms.com/help/mvpstats/distributions/SkewnessCriticalValues
#
#     with open(os.path.join(folder,'kurtosis.txt'), 'wt') as f:
#         printListOfListPerIteration(makeCharactOutOfEndDifferencesList(a,STAT_KURTOSIS), f) #critical low -0.54	high 0.79  http://mvpprograms.com/help/mvpstats/distributions/KurtosisCriticalValues
#

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







#
if __name__ == '__main__':
    test()



