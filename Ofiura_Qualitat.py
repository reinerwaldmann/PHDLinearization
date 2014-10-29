__author__ = 'vasilev_is'

import math

import numpy as np

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
        diflist.append(np.array(measpoint['y'])-func(measpoint['x'],b,c))
    return np.average(diflist), np.var(diflist), math.sqrt(np.var(diflist))

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


