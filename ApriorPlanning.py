__author__ = 'vasilev_is'

import copy
import random
import math

import numpy as np


def replaceInList (list, i, x):
    """
    :param list: список
    :param i: индекс элемента
    :param x: переменная, на которую заменяем элемент
    :return: Возвращает список list, в котором i-й элемент заменяется на x
    """
    l=copy.deepcopy(list)
    l[i]=x
    return l

def uniformVector(xstart, xend):
    """
    :param xstart: вектор начала диапазона
    :param xend: вектор конца диапазона
    :return: вектор, являющийся частью равномерного распределения между xstart и xend
    """
    res= [0 for i in range(len(xstart))]
    for i in range(0, len(xstart)):
        res[i]= random.uniform (xstart[i], xend[i])
    return res
    #Можно сделать map для списка списков сразу (или просто для нескольких списков) на основе этой функции



def makeUniformExpPlan(xstart:list, xend:list, N:int):
    """
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: равномерный план эксперимента
    """

    res=list()

    xstartnp=np.array(xstart)
    xendnp=np.array(xend)

    xstep = (xendnp-xstartnp)/N

    for i in range(N):
        res.append({"x":list(xstart+i*xstep) })

    return res


def makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    """
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор
    """
    res=list()
    for i in range(0, N):
        res.append({"x": uniformVector(xstart, xend)})
    return res



def makeMeasAccToPlan(func, expplan:list, b:list, c:dict, ydisps:list=None, n=1, outfilename="", listOfOutvars=None):
    """

    :param func: векторная функция
    :param expplan: план эксперимента
    :param b: вектор b
    :param c: вектор c
    :param ydisps: ковариационная матрица y
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """

    res = copy.deepcopy(expplan)

    for i in range(len(res)):
        y=func(res[i]["x"],b,c)
        #Внесём возмущения:
        if not ydisps==None:
            for i in range (len(y)):
                y[i]=random.normalvariate(y[i], math.sqrt(ydisps[i]))

        res[i]["y"]=y
    return res







def test():
    xstart=[0,20,100]
    xend=[10,30,200]
    N=10

    print (makeUniformExpPlan(xstart, xend, N))



test()