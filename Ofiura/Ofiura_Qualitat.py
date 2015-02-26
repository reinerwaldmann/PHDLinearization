__author__ = 'vasilev_is'
import copy
import math
import os

import numpy as np
from prettytable import PrettyTable
"""
logTruthness (measdata:list, b:list, Ve,  func, c):
    Считает логарифм функции правдоподобия для известной ковариационной матрицы Ve - ошибок экспериментальных данных
    :param measdata: измеренные данные
    :param b: оценка вектора коэффициентов
    :param Ve: ковариационная матрица ошибок экспериментальных данных
    :param func: callable функция x,b,c, возвращает значение y
    :param c: словарь дополнительных переменных
    :return: среднее значение логарифма функции правдоподобия, дисперсию по выборке экспериментальных данных, стандартное отклонение

def averageDif(measdata:list, b:list, Ve,  func, c):
    :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param b: вектор коэфф
    :param Ve:  Ve ковар. матрица измеренных данных
    :param funcf callable функция, параметры по формату x,b,c
    :param c словарь дополнительных постоянных
    :return: среднее, дисперсия, стандартное отклонение

def getQualitat(measdata:list, b:list, Ve,  func, c):
    :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param b: вектор коэфф
    :param Ve:  Ve ковар. матрица измеренных данных
    :param funcf callable функция, параметры по формату x,b,c
    :param c словарь дополнительных постоянных
    :return: Среднее логарифма правдоподобия Дисперсия лп Сигма лп Среднее остатков Дисп. остатков Сигма остатков

"""

foldername=None


def convertToQualitatStandart (gknux, funcf, jacf,  measdata:list, c, Ve, name=''):
    """
    Принимая на вход gknux по старому стандарту (tuple) и некоторые другие параметры, возвращает словарь с gknux - параметрами
    и показателями качества, то есть новый стандарт

    Спецификация:
    Тип данных - словарь
    Поля и описание -
     'b'            оцененный вектор коэффициентов
     'numiter'      число итераций
     'log'          лог
     'Sklist'       список значений объектной функции
     'Sk'           последнее значение объектной функции
     'AvLogTruth'   среднее по точкам эксперимента значение логарифма правдоподобия
     'DispLT'       дисперсия логарифма правдоподобия
     'SigmaLT'      среднеквадратическое отклонение логарифма правдоподобия
     'AvDif'        средний остатков
     'DispDif'      дисперсия остатков
     'SigmaDif'     СКВ остатков
     'Diflist'      список остатков
     'name'         название метода
     'Vb'           ковариационная матрица оценённого вектора b
     'VbSigmas'     сигмы из ковариационной матрицы оценённого вектора b (корни квадратные из диагонали)

    :param gknux- результат оценочной функции по старому стандарту
    :param funcf- указатель на функцию
    :param jacf - указатель на функцию, возвращ. якобиан
    :param measdata- измеренные данные
    :param c
    :param Ve
    :param name - название метода


    """
    names = [' b', 'numiter', 'log' , 'Sklist', 'Sk']
    gknuxdict = dict(zip (names, list(gknux)))
    #    names=['AvLogTruth','DispLT', 'SigmaLT', 'AvDif', 'DispDif', 'SigmaDif', 'Diflist']
    rs= dict(list(getQualitatDict(measdata, gknuxdict[' b'], Ve,  funcf, c).items()) + list(gknuxdict.items()))
    rs['name']=name
    Vb = countVbForMeasdata(gknux[0],  c, Ve, jacf, measdata)
    rs['Vb'] = Vb
    sigmas=list()
    for i in range (Vb.shape[0]):
        sigmas.append(math.sqrt(Vb[i][i]))
    rs['VbSigmas'] = sigmas
    return rs


def printQualitatStandart(gknux_dict:dict):
    """
    Выводит таблицу показателей качества оценки и логи
    USES PRETTYTABLE
    :param gknux_dict результаты оценки с показателями качества по стандарту 2
    """
    g=gknux_dict
    tg=copy.copy(g) #для вывода в таблицу

    del tg['log']
    del tg['Sklist']
    del tg['Diflist']
    del tg['name']
    del tg['Vb']



    klist = sorted(list(tg.keys()))
    vallist = [tg[i] for i in klist]


    num=4
    klists=[klist[i::num] for i in range(num)]
    vallists=[vallist[i::num] for i in range(num)]


    t=PrettyTable (klists[0])
    t.add_row(vallists[0])

    t2=PrettyTable (klists[1])
    t2.add_row(vallists[1])


    print ("\n Results of method {0} \n".format(g['name']))
    print (t)
    print (t2)
    print ("Vb:\n {0}  ".format (g['Vb']))
    print ("\nLog messages: {0} \n ".format (g['log']))


def printQualitatNeat(measdata:list, b:list, Ve,  func, c, jac):
    """
    Выводит таблицу показателей качества оценки
    USES PRETTYTABLE
    :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param b: вектор коэфф
    :param Ve:  Ve ковар. матрица измеренных данных
    :param funcf callable функция, параметры по формату x,b,c
    :param c словарь дополнительных постоянных
    :return: Среднее логарифма правдоподобия Дисперсия лп Сигма лп Среднее остатков Дисп. остатков Сигма остатков
    """

    t=PrettyTable (['Среднее логарифма правдоподобия','Дисперсия лп', 'Сигма лп', 'Среднее остатков', 'Дисп. остатков', 'Сигма остатков'])
    t.add_row(list(logTruthness (measdata, b, Ve,  func, c))+list(averageDif(measdata, b, Ve,  func, c))[:-1:]  )
    print('Показатели качества оценки')
    print (t)

    print('Матрица Vb')

    Vb=countVbForMeasdata(b,  c, Ve, jac, measdata)
    print (Vb)

    print('Parameter Sigmas')
    for i in range (Vb.shape[0]):
        print (math.sqrt(Vb[i][i]))





def logTruthness (measdata:list, b:list, Ve,  func, c):
    """
    Считает логарифм функции правдоподобия для известной ковариационной матрицы Ve - ошибок экспериментальных данных
    :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param b: оценка вектора коэффициентов
    :param Ve: ковариационная матрица ошибок экспериментальных данных
    :param func: callable функция x,b,c, возвращает значение y
    :param c: словарь дополнительных переменных
    :return: среднее значение логарифма функции правдоподобия, дисперсию по выборке экспериментальных данных, стандартное отклонение
    """
    S=list()

    for i in range(len(measdata)):
        measpoint = measdata[i]
        dif=np.array(measpoint['y'])-np.array(func(measpoint['x'],b,c))

        S.append(np.dot(np.dot(dif.T, np.linalg.inv(Ve)), dif))


    K=Ve.shape[0] #число откликов
    M=len(b) #число коэффициентов
    N=len(measdata)
    shift=K*N/(K*N-M)

    Average=np.average(S)*shift
    Disp = np.var(S)*shift*shift

    return Average, Disp, math.sqrt(Disp)

def averageDif(measdata:list, b:list, Ve,  func, c):
    """
    :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param b: вектор коэфф
    :param Ve:  Ve ковар. матрица измеренных данных
    :param funcf callable функция, параметры по формату x,b,c
    :param c словарь дополнительных постоянных
    :return: среднее, дисперсия, стандартное отклонение
    """

    diflist=list()
    diflistNoAbs=list()


    for measpoint in measdata:
        dif=np.array(measpoint['y'])-func(measpoint['x'],b,c)

        diflist.append(np.abs(dif ))
        diflistNoAbs.append(dif)

    return np.average(diflist), np.var(diflist), math.sqrt(np.var(diflist)), diflistNoAbs

def getQualitat(measdata:list, b:list, Ve,  func, c):
    """
    :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param b: вектор коэфф
    :param Ve:  Ve ковар. матрица измеренных данных
    :param funcf callable функция, параметры по формату x,b,c
    :param c словарь дополнительных постоянных
    :return: Среднее логарифма правдоподобия Дисперсия лп Сигма лп Среднее остатков Дисп. остатков Сигма остатков
    """
    return "Среднее логарифма правдоподобия Дисперсия лп Сигма лп Среднее остатков Дисп. остатков Сигма остатков\n", logTruthness (measdata, b, Ve,  func, c), averageDif(measdata, b, Ve,  func, c)


def getQualitatDict(measdata:list, b:list, Ve,  func, c):
    names=['AvLogTruth','DispLT', 'SigmaLT', 'AvDif', 'DispDif', 'SigmaDif', 'Diflist']
    values = list(logTruthness (measdata, b, Ve,  func, c))+list(averageDif(measdata, b, Ve,  func, c))
    return dict(zip (names, values))





def countVbForMeasdata(b:list,  c:dict, Ve, jac, measdata):
    """
    Для неявных функций актуальным является значение y. В некоторых случаях (напр. при последовательном планировании)
    y уже просчитан, и нет нужды его считать снова, задавая функцию.
    :param expplan: план эксперимента
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :param measdata: данные измерений
    :return: значение определителя для данного плана эксперимента
    """
    G=np.zeros((len(b),len(b))) #матрица G

    for point in measdata:
        jj=jac(point['x'], b, c, point['y'])
        #G+=jj*np.linalg.inv(Ve)*jj.T
        G+=np.dot ( np.dot(jj.T, np.linalg.inv(Ve)), jj)

        #G+=np.dot(jj.T, jj)
    try:
        return np.linalg.inv(G)
    except BaseException as e:
        print('Fatal error in countVbForMeasdata: ',e)
        print('b vector=',b)
        print('current point=',point)
        print('G=',G)
        exit(0)





def printGKNUNeat(gknu):
    """
    Выводит таблицу результатов работы метода оценки
     USES PRETTYTABLE
     """
    t=PrettyTable(['b','Количество итераций', 'log', 'Skmu'])

    if type(gknu)==tuple:
        t.add_row([gknu[0],gknu[1],gknu[2],gknu[4]])
    elif type(gknu)==dict:
        t.add_row([gknu['b'],gknu['numiter'],gknu['log'],gknu['Sk']])
    else:
        print  ('analyseDifList erroneous arg')
        return


    print('Данные оценки')
    print (t)


def printSeqPlanData(seq):
    """
    Выводит таблицу результатов работы метода последовательного планирования
     USES PRETTYTABLE
     """
    t=PrettyTable(['Число итераций - число добавленных точек','Sk', 'Логи', 'detVb'])
    t.add_row([seq[1], seq[2], seq[5], seq[8]])
    print('Данные работы последовательного плана')
    print (t)



def analyseDifList(arg, plot=True, imagename='img' ): #числовой режим, вектора не поддерживаются
    """
    Функция, выводит гистограмму остатков
    в будущем будет возвращать True, если распределение нормальное и false, если таки нет
    пока возвращает среднее
    :param arg
    :param plot
    """
    global foldername


    if type(arg)==dict:
        diflist = arg['Diflist']
    elif type(arg)==list:
        diflist = arg
    else:
        print  ('analyseDifList erroneous arg')
        return None

    superlist = list(map(lambda x: x[0], diflist)) #костыль для тех случаев, когда на вход всё одно вектор, но единичной длины

    #1. Определение нормальности распределения
    #пока не имплементим

    # учитываем, что всё это какбе с векторами
    # с векторами можно работать как - 1. рассматривая каждый компонент по отдельности


    if plot:
        #посмотреть диаграмму остатков у лучшей оценки
        import matplotlib.pyplot as plt
        plt.hist(superlist, 25, label=imagename)

        if foldername:
            try:
                os.makedirs(foldername)
            except OSError:
                if not os.path.isdir(foldername):
                    raise
            plt.savefig(foldername+imagename+'.png')
        else:
            plt.show()

    return np.average(superlist), True











