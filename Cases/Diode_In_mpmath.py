from mpmath.calculus.optimization import MDNewton

__author__ = 'vasilev_is'

import math

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import mpmath as mpm


#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit
#Опция: MPMATH USED


#Таблица переменных:

#Установка точности для mpmath. В проекте принято 50
mpm.mp.dps=50
mpm.pretty = True


#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=mpm.mpf('0.02586419') #подогнанное по pspice


#Diode_In_mpmath

def func_Diode_In_mpmath(y,x,b,c):
    """
    [Реестровая] +
    Уравнение Кирхгофа
    :param y:
    :param x:
    :param b:
    :param c:
    :return:
    """
    global FT
    mm=b[0]*(mpm.exp((x[0]-y[0]*b[2])/(FT*b[1]))-1)-y[0]
    return [mm]



def solver_Diode_In_mpmath (x,b,c=None):
    """
    [Реестровая] +
    :param x:
    :param b:
    :param c:
    :return:
    """
    global numnone
    global FT
    #первая часть решения: находим корень с небольшой точностью через numpy
    dfdy=lambda y,x,b,c=None: np.array ([[ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]])
    func = lambda y,x,b,c: [np.float(x) for x in func_Diode_In_mpmath(y,x,b,c)]
    solvinit=[0.5]
    try:
        #solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        solx=optimize.root(func, solvinit, args=(x,b,c), method='lm').x
        #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        #print ('diodeResistorIMPLICITfunction: Error in findroot=',e)
        numnone+=1
        print ("solver: ERROR: "+e.__str__())
        return [None]

    #вторая часть решения: находим корень с увеличенной  точностью через mpmath
    funcMPM = lambda y:  func_Diode_In_mpmath ([y],x,b,c)
    dfdyMPM=lambda y: mpm.matrix ([[ -1 - b[0]*b[2]*mpm.exp((-b[2]*[y][0] + x[0])/(FT*b[1]))/(FT*b[1])]])
    solvinitMPM=solx[0]
    try:
        precsolx=mpm.calculus.optimization.findroot(mpm.mp, f=funcMPM, x0=solvinitMPM, solver=MDNewton, multidimensional=True, J=dfdyMPM)
    except BaseException as e:
        numnone+=1
        print ("solver MPM: ERROR: "+e.__str__())
        return [None]

    return precsolx






def test_Diode_In_mpmath():
    b=[mpm.mpf('1.238e-14'), mpm.mpf('1.3'), mpm.mpf('1')]
    #0.0026638081177255196
    rng=mpm.arange('0.01','2','0.01')
    #снимем ВАХ
    resrng=[solver_Diode_In_mpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
#    resrngorig=[casesDiode.diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    b=[mpm.mpf('1.238e-14'), mpm.mpf('1.3'), mpm.mpf('100')]
    resrng1=[solver_Diode_In_mpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
 #   plt.plot(rng , resrngorig, label='r=0')
    plt.plot(rng , resrng, label='r=1000')
    plt.plot(rng , resrng1, label='r=3000')
    plt.legend(loc='upper left')
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()


test_Diode_In_mpmath()

