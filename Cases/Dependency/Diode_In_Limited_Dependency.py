
__author__ = 'reiner'

import math

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

import Ofiura.Ofiura_Qualitat as o_q
import Ofiura.Ofiura_EstimationLimited as o_el
import Ofiura.Ofiura_Estimation  as o_e
import Ofiura.Ofiura_ApriorPlanning as o_ap
import Ofiura.Ofiura_planning as o_p



#Этот файл полностью повторяет по сути Diode_In_Limited, но адаптирован к эксперименту Dependency


#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit

#Таблица переменных:


#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=0.02586419 #подогнанное по pspice
foldername='cachedPlans'



"""
Экстрагируем три параметра диода: N (коэфф. неидеальности), Is, R(омическое сопротивление, включенное последовательно с диодом)
"""

#Внимание - надо исключать из плана эксперимента неверно оценённые точки. Неверными признаются те, у которых результат совпадает с init
#В этом случае рабочая функция должна возвращать none, а создатель плана - вежливо вытряхивать точку из measdata


def func_Diode_In_Limited(y,x,b,c=None):
    FT=0.02586419 #подогнанное по pspice
    mm=float(b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0])
    return [mm]




def solver_Diode_In_Limited (x,b,c=None):
    """
    Возвращает ток, текущей по цепи

    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, последовательное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода

    ===ОГРАНИЧЕНИЯ===
    #Сочетание боьльших напряжений  (>1В) и нулевого сопротивления даёт идиотский результат (единицу, то есть метод не отрабатывает)
    Аналогичное нехорошее поведение может повстречаться и при прочих сочетаниях. Этот момент должно иметь в виду

    :return:s
    """

    global numnone
    global resultsOfEstimation

    # V=x[0] #напряжение на диоде
    # Is=b[0]
    # N=b[1]
    # R=b[2]
    global FT

    dfdy=lambda y,x,b,c=None: np.array ([[ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]])

    solvinit=[1]

    try:
        solx=optimize.root(func_Diode_In_Limited, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        print ('diodeResistorIMPLICITfunction: Error in findroot=',e)
        numnone+=1
        return None

    if solx-solvinit==[0]*len(solx):
         numnone+=1
         return None

    return solx


def jac_Diode_In_Limited (x,b,c,y):
    global FT


    dfdb=lambda y,x,b,c: np.matrix( [ [math.exp((-b[2]*y[0] + x[0])/(FT*b[1])) - 1,
                                      -b[0]*(-b[2]*y[0] + x[0])*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                      -b[0]*y[0]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]     ])


    dfdy=lambda y,x,b,c: np.matrix ([ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])])

    jacf=np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


    return jacf



def plotPlanAndMeas2D(measdata):

    planplot1=[x['x'][0] for x in measdata]
    measplot1=[x['y'][0] for x in measdata]
    plt.plot(planplot1, measplot1,  'ro')
    plt.ylabel('value')
    plt.xlabel('point')
    plt.grid()
    plt.show()




def extraction_Diode_In_Limited():
    """
    пробуем экстрагировать коэффициенты из модели диода
    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, параллельное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода
    :return:
    """

    global FT
    global foldername

    #возвращает значение y
    funcf=solver_Diode_In_Limited
    jacf = jac_Diode_In_Limited
    #теперь попробуем сделать эксперимент.
    c={}
    Ve=np.array([ [0.00000001] ]  )
    btrue=[1.238e-14, 1.8, 100]
    #btrue=[1.5e-14, 1.75, 150]

    bstart=np.array(btrue)-np.array(btrue)*0.3
    bend=np.array(btrue)+np.array(btrue)*0.3002
    #binit=np.array(btrue)-np.array(btrue)*0.1

    print('conditions:')
    print(bstart)
    print(bend)


    binit=[1.1e-14, 1.5, 90]
    xstart=[0.001]
    xend=[2]
    N=20
    print("performing aprior plan:")

#примитивная попытка автоматизировать, риальни надо кешировать в файл под хешем параметров

    import os
    filename =foldername+'/'+os.path.basename(__file__).replace('.py','_plan')
    try:
        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:
        oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbose=True)[1]
        o_p.writePlanToFile(oplan, filename)

#    получаем измеренные данные
    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve )
#     #чертим эти данные
#     #o_pl.plotPlanAndMeas2D(measdata, 'Aprior Disp{0} measdata'.format(Ve))
#
#     #оценка


    #grandCountGN_UltraX1_Limited_wrapper (funcf, jacf,  measdata:list, binit:list, bstart:list, bend:list, c, A, NSIG=50, NSIGGENERAL=50, implicit=False, verbose=False, verbose_wrapper=False):

    gknuxlim = o_el.grandCountGN_UltraX1_Limited_wrapper(funcf,jacf,measdata,binit,bstart,bend, c, implicit=True, verbose=False, verbose_wrapper=False )
    gknux = o_e.grandCountGN_UltraX1(funcf, jacf, measdata, binit, c, implicit=True)

    gknuxlim2=o_q.convertToQualitatStandart (gknuxlim, funcf, jacf,  measdata, c, Ve, name='Limited Count Aprior')
    gknux2=o_q.convertToQualitatStandart (gknux, funcf, jacf,  measdata, c, Ve, name='Normal Count Aprior')

    o_q.printQualitatStandart (gknuxlim2)
    o_q.printQualitatStandart (gknux2)






#extraction_Diode_Irev_Limited()