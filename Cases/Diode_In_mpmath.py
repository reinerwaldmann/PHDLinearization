from mpmath.calculus.optimization import MDNewton

__author__ = 'vasilev_is'

import math

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import mpmath as mpm
import Ofiura.Ofiura_planningMpmath as o_pmpm
import Ofiura.Ofiura_EstimationMpmath as o_empm

#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit
#Опция: MPMATH USED

#Таблица переменных:
#Установка точности для mpmath. В проекте принято 50
# mpm.mp.dps=40
# mpm.pretty = True

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
        numnone+=1
        print ("solver: ERROR: "+e.__str__())
        return [None]

    #вторая часть решения: находим корень с увеличенной  точностью через mpmath
    funcMPM = lambda y:  func_Diode_In_mpmath ([y],x,b,c)
    dfdyMPM=lambda y: mpm.matrix ([[ -1 - b[0]*b[2]*mpm.exp((-b[2]*[y][0] + x[0])/(FT*b[1]))/(FT*b[1])]])
    solvinitMPM=solx[0]
    try:
        precsolx=mpm.calculus.optimization.findroot(mpm.mp, f=funcMPM, x0=solvinitMPM, solver=MDNewton, multidimensional=True, J=dfdyMPM, verify=False)
    except BaseException as e:
        numnone+=1
        print ("solver MPM: ERROR: "+e.__str__())
        return [None]
    return precsolx


def jac_Diode_In_mpmath (x,b,c,y):
    """
    [Реестровая]
    непосредственная проверка невозможна


    :param x:
    :param b:
    :param c:
    :param y:
    :return:
    """
    global FT
    dfdbMPM=lambda y,x,b,c: mpm.matrix( [ [mpm.exp((-b[2]*y[0] + x[0])/(FT*b[1])) - 1,
                                          -b[0]*(-b[2]*y[0] + x[0])*mpm.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                          -b[0]*y[0]*mpm.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]     ])
    dfdyMPM=lambda y,x,b,c: mpm.matrix ([[ -1 - b[0]*b[2]*mpm.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]])
    jacf=dfdyMPM(y,x,b,c)**-1 * dfdbMPM(y,x,b,c)
    return jacf


def test_Diode_In_mpmath():
    b=[mpm.mpf('1.238e-14'), mpm.mpf('1.3'), mpm.mpf('1')]
    c=None
    #0.0026638081177255196
    rng=mpm.arange('0.01','2','0.01')
    xstart=[mpm.mpf('0.00001')]
    xend=[mpm.mpf('2')]
    Ve=np.array([ [0.0001] ]  )

    #снимем ВАХ
#    resrng=[solver_Diode_In_mpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

    #resrng=o_pmpm.makeMeasOneDot_lognorm(solver_Diode_In_mpmath, x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    resrng=[o_pmpm.makeMeasOneDot_lognorm(solver_Diode_In_mpmath, [x],b,c,Ve) for x in rng]

    resrng1=[solver_Diode_In_mpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
 #   plt.plot(rng , resrngorig, label='r=0')
    plt.plot(rng , resrng, label='r=1000')
    plt.plot(rng , resrng1, label='r=3000')
    plt.legend(loc='upper left')
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()



def extraction_Diode_In_mpmath():
    """
    пробуем экстрагировать коэффициенты из модели диода
    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, параллельное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода
    :return:
    """
    btrue=[mpm.mpf('1.238e-14'), mpm.mpf('1.3'), mpm.mpf('10')]
    c=None
    global FT
    funcf=solver_Diode_In_mpmath
    jacf = jac_Diode_In_mpmath
    #теперь попробуем сделать эксперимент.
    Ve=np.array([ [1e-7] ]  )
    bstart=[mpm.mpf('1.0e-14'), mpm.mpf('1.0'), mpm.mpf('9')]
    bend=[mpm.mpf('1.5e-14'), mpm.mpf('1.5'), mpm.mpf('14')]
    #binit=[mpm.mpf('1.1e-14'), mpm.mpf('1.1'), mpm.mpf('11')]

    #binit=mpm.matrix([['1.1e-14', '1.1', '11.1']]).T
    binit=btrue

    xstart=[mpm.mpf('0.0001')]
    xend=[mpm.mpf('1.5')]

    N=100
    NAprior=20

    unifplan = o_pmpm.makeUniformExpPlan(xstart, xend, N)
    unifmeasdata = o_pmpm.makeMeasAccToPlan_lognorm(funcf, unifplan, btrue, c, Ve)
    print (unifmeasdata)
    #print (unifmeasdata[0]['y'][0])

    gknux = o_empm.grandCountGN_UltraX1_mpmath(funcf, jacf, unifmeasdata, binit,c, NSIG=100,implicit=True, verbose=True)
    print (gknux)



mpm.mp.dps=5
extraction_Diode_In_mpmath()







#test_Diode_In_mpmath()





