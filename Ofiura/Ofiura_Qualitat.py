__author__ = 'vasilev_is'

import math

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
    for measpoint in measdata:
        diflist.append(np.abs(np.array(measpoint['y'])-func(measpoint['x'],b,c)))
    return np.average(diflist), np.var(diflist), math.sqrt(np.var(diflist)), diflist

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









def printQualitatNeat(measdata:list, b:list, Ve,  func, c):
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


def printGKNUNeat(gknu):
    """
    Выводит таблицу результатов работы метода оценки
     USES PRETTYTABLE
     """
    t=PrettyTable(['b','Количество итераций', 'log', 'Skmu'])
    t.add_row([gknu[0],gknu[1],gknu[2],gknu[4]])
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




