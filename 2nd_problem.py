

__author__ = 'reiner'
"""
В файле реализованы функции, потребные для решения проблем во втором задании.
Краткое описание:

1. Генерируем значение r1,r2  +
2. Применяем к ним формулу делителя  +
3. Вырезаем нужную область +
4. Перенсоим вырезанный диапазон на область r1-r2  +
5. В этой области строим прямоугольник

"""

import copy

import numpy as np
import matplotlib.pyplot as plt

import sequence_generation as sg


def filterfunc (input:dict, diaps:dict):
    """
    На входе - словарь входных значений, словарь их диапазонов
    На выходе - подошёл (1) или не подошёл (0)

    @input - словарь вх знач
    @diaps - словарь диапазонов
    """

    # {"a": [1, 0] , "b":[10, 1000000]}

    for key, diap in diaps.items():
        if not diap[0]<input[key]<diap[1]:
            return 0
    return 1


def filterseq (inputlist:list, diaps:dict):
    """
    Фильтрует и определяет процент годных
    На входе:
    @inputlist - список словарей входных значений
    @diaps - словарь диапазонов

    На выходе:
    Фильтрованная последовательность, процент выхода годных (процент, который составляет выходная последовательность от входной)

    """

    filterfuncl=lambda x: filterfunc(x, diaps)
    filtered=list(filter(filterfuncl, inputlist))

    return filtered, 100*len(filtered)/len(inputlist)    # В результат попадают только те элементы x, для которых x < 5 истинно


def showonlyoutputseq(inputseq, var):
    """
    Получает последовательность словарей,
    возвращает последовательность значений лишь одной переменной

    @niputseq - вх после
    @var - переменная, значения которой мы извлекаемs
    """
    #отображаем только ksu
    filterksu = lambda x: x[var]
    seq=copy.deepcopy(inputseq)

    return list (map (filterksu, seq))

def getellipse2D (seq, strx, stry):
    """
    1. Получает все параметры эллипса. А именно:
     матожидание обоих параметров, дисперсии обоих параметров, коэффициент корреляции обоих параметров.
     этими параметрами определяется уравнение эллипса

    """
    xseq=showonlyoutputseq(seq, strx)
    yseq=showonlyoutputseq(seq, stry)

    mx=np.mean(xseq)
    my=np.mean(yseq)

    dx=np.var(xseq)
    dy=np.var(yseq)

    sigmax=np.sqrt(dx)
    sigmay=np.sqrt(dy)


    r=np.corrcoef(xseq, yseq)[0][1]






    # print (r)
    #variance - это дисперсия, то есть сигма в квадрате








#test area

funcstrdict= {"ksu2":"(r2)/(r1+r2)"}


xvectorlistsdict = {"r1":[20], "r2":[30]}

spreadvarslist  = ["r1", "r2"]


V1=np.array                   ( [[4, 2],
                                [3, 6]])

V=np.array                   ( [[4, 0],
                                [0, 6]])

resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 10000)



diaps={"ksu2":[0.6, 0.61]}
filteredresdict=filterseq (resdict, diaps)[0]

seqr1=showonlyoutputseq (resdict, "r1")
seqr2=showonlyoutputseq (resdict, "r2")

seqr1f=showonlyoutputseq (filteredresdict, "r1")
seqr2f=showonlyoutputseq (filteredresdict, "r2")


getellipse2D (filteredresdict, "r1", "r2")


exit(0)

plt.figure(1) # Here's the part I need, but numbering starts at 1!
plt.xlabel('seqr1')
plt.ylabel('seqr2')
#plt.axis([0.94, 0.97, 0.86, 0.92])
plt.plot(seqr1, seqr2, 'ro')

plt.figure(2) # Here's the part I need, but numbering starts at 1!
plt.xlabel('seqr1')
plt.ylabel('seqr2')

plt.plot(seqr1f, seqr2f, 'ro')




plt.show()


# print ("Исходная последовательность ksu2")
# print (showonlyoutputseq(resdict, "ksu2"))
# print (resdict)
# print ("Фильтрованная последовательность ksu2")
# print (showonlyoutputseq(filteredresdict, "ksu2"))
# print (filteredresdict)
#
#
