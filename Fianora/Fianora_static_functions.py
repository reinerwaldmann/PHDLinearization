
__author__ = 'vasilev_is'

import copy
import random
import sys

"""
replaceInList (list, i, x):
    :param list: список
    :param i: индекс элемента
    :param x: переменная, на которую заменяем элемент
    :return: Возвращает список list, в котором i-й элемент заменяется на x
 uniformVector(xstart, xend):
    :param xstart: вектор начала диапазона
    :param xend: вектор конца диапазона
    :return: вектор, являющийся частью равномерного распределения между xstart и xend

doublesearch (xstart, xend, xinit, function):
    Реализует метод двойного поиска с ограничениями. Внимание: функция требует правки в части штрафов в зависимости от значений, которые может принимает объектная функция.
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param xinit:начальное значение x (вектор)
    :param function: объектная функция, принимает на вход x
    :return: оптимизированное значение x
"""




class EstimationContext():
    def __init__(self, bstart, bend, btrue, binit, xstart, xend, Ve, N):
        self.bstart = bstart
        self.bend = bend
        self.btrue = btrue
        self.binit = binit
        self.xstart = xstart
        self.xend = xend
        self.Ve = Ve
        self.N = N




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

def appendToList (list, x):
    """
    :param list: список
    :param x: переменная, которую добавляем
    :return: Возвращает список list, в который добавляется x
    """
    l=copy.deepcopy(list)
    l.append(x)
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
    #Можно сделать map для списка списков сразу (или просто для нескольких списков) на основе этой функции  #


def rangomNormalvariateVector (xstart, xend):

    res = [0 for i in range(len(xstart))]
    middle = [0 for i in range(len(xstart))]
    sigma = [0 for i in range(len(xstart))]

    for i in range(0, len(xstart)):
        middle[i]=xstart[i]+(xend[i]-xstart[i])/2
        sigma[i]=(xend[i]-middle[i])/3
        res[i]= random.normalvariate (middle[i], sigma[i])

    return res


def doublesearch (xstart, xend, xinit, function):
    """
    Метод двойного поиска с ограничениями
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param xinit:начальное значение x (вектор)
    :param function: объектная функция, принимает на вход x
    :return: оптимизированное значение x
    """
    x=copy.deepcopy(xinit)
    #xcurr=copy.deepcopy(x)
    for i in range (len(x)):
        step=(xend[i]-xstart[i])/3 #задание начального шага
        brval=10000000 #критичное число итераций, после которого надо выйти из цикла и сообщить об ошибке
        while(1):
            xcurr=copy.deepcopy(x)
            d3=function(xcurr)
            xcurr[i]=x[i]+step  #пробуем шаг вперёд
            d4=function(xcurr) if xcurr[i]<xend[i] else sys.maxsize #штраф при выходе из диапазона
            xcurr[i]=x[i]-step #пробуем шаг назад
            d5=function(xcurr) if xcurr[i]>xstart[i] else sys.maxsize #штраф при выходе из диапазона

            # print ('head')
            # print (d3,d4,d5, step, x, i)
            # print ('tail')


            if (d3-d4)*(d3-d5)>0: #
                step=step/2       #уменьшение шага вполовину
                if step<(xend[i]-xstart[i])/1000: #если шаг меньше 1/1000, значит этот член вектора прооптимизировали, переходим к следующему
                    break
            else:
                if (d4>d5):
                    x[i]=x[i]-step #шаг назад
                else:
                    x[i]=x[i]+step #шаг вперёд

            brval-=1
            if not brval:
                print ("doublesearch: exit due to max. number of iteration archieved, possible ERROR")
                break

    return x


import time

class Profiler(object):
    def __enter__(self):
        self._startTime = time.time()

    def __exit__(self, type, value, traceback):
        print ("Elapsed time: {:.3f} sec".format(time.time() - self._startTime))