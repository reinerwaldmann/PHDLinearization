__author__ = 'vasilev_is'

import math
import os
import decimal as dc

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import Ofiura.Ofiura_Estimation as o_e
import Ofiura.Ofiura_ApriorPlanning as o_ap
import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Qualitat as o_q




#Задаём имя папки, куда складывать нарисованные изображения
genfoldername = '/home/reiner/results_estimation/'
subfoldername = os.path.basename(__file__).replace('.py','_results/')
o_q.foldername = genfoldername + subfoldername



#Part1: прямая ветвь ВАХ диода
#В файле предполагается реализовать 2 подхода:
# - переход к decimal и т.е. увеличение точности
# - переход к задаче с ограничениями


#Внимание: здесь используется не преобразованный вариант с резистором  в цепи




#Стандарт: implicit
#Таблица переменных:

#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=np.float128(0.02586419) #подогнанное по pspice


def func_DiodeV4_partIn_Decimal(y,x,b,c=None):
    global FT


    mm=np.float128(b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0])
    return [mm]




def solver_DiodeV4_partIn_Decimal (x,b,c=None):
    """
    Возвращает ток, текущей по цепи

    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, последовательное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода

    ===ОГРАНИЧЕНИЯ===
    #Сочетание боьльших напряжений  (>1В) и нулевого сопротивления даёт идиотский результат (единицу, то есть метод не отрабатывает)
    Аналогичное нехорошее поведение может повстречаться и при прочих сочетаниях. Этот момент должно иметь в виду

    :return:
    """

    global numnone
    global resultsOfEstimation
    global FT

    # V=x[0] #напряжение на диоде
    # Is=b[0]
    # N=b[1]
    # R=b[2]


    dfdy=lambda y,x,b,c=None: np.array ([[dc.Decimal( -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]))]])




    #function=lambda y,x,b,c=None: [b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0]] #в подготовленном для взятия производной виде



    solvinit=[dc.Decimal(1)]

    try:
        solx_source=optimize.root(func_DiodeV4_partIn_Decimal, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        solx = [dc.Decimal(x) for x in solx_source]
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





def jac_DiodeV4_partIn_Decimal (x,b,c,y):
    FT=0.02586419 #подогнанное по pspice


    dfdb=lambda y,x,b,c: np.matrix( [ [math.exp((-b[2]*y[0] + x[0])/(FT*b[1])) - 1,
                                      -b[0]*(-b[2]*y[0] + x[0])*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                      -b[0]*y[0]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]     ])


    dfdy=lambda y,x,b,c: np.matrix ([ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])])

    #возвращает структурную матрицу
    #jacf=lambda x,b,c,y: jjacf(x,b,c,y,dfdb,dfdy)


    jacf=np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


    return jacf



def test__DiodeV4_partIn_Decimal():

    b=[1.238e-14, 1.3, 0    ]
    #0.0026638081177255196
    rng=np.arange(0.01,2,0.01)
    #снимем ВАХ
    resrng=[solver_DiodeV4_partIn_Decimal([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
#    resrngorig=[casesDiode.diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    b=[1.238e-14, 1.8, 3000]
    resrng1=[solver_DiodeV4_partIn_Decimal ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
 #   plt.plot(rng , resrngorig, label='r=0')
    plt.plot(rng , resrng, label='r=1000')
    plt.plot(rng , resrng1, label='r=3000')
    plt.legend(loc='upper left')
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()




def extraction_DiodeV4_partIn_Decimal(plot=True):
    """
    пробуем экстрагировать коэффициенты из модели диода
    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, параллельное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода
    :return:
    """

    #возвращает значение y
    funcf=solver_DiodeV4_partIn_Decimal
    jacf = jac_DiodeV4_partIn_Decimal
    #теперь попробуем сделать эксперимент.
    c={}
    Ve=np.array([ [0.00000001] ]  )
    btrue=[1.238e-14, 1.8, 100]
    #btrue=[1.5e-14, 1.75, 150]

    bstart=np.array(btrue)-np.array(btrue)*0.3
    bend=np.array(btrue)+np.array(btrue)*0.3
    binit=[1.0e-14, 1.7, 70]
    #binit=[1.238e-14, 1.8, 100]

    xstart=[0.001]
    #xend=[20,60]
    xend=[1.5]
    N=100 #для неаприорных планов
    NAprior=20 #для априорных планов


     #Получаем априорный план
    print("performing aprior plan:")
    #блок кеширования априорного плана в файл
    filename ='cachedPlans/'+os.path.basename(__file__).replace('.py','_plan')
    try:
        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:

        oplan=o_ap.grandApriornPlanning (xstart, xend, NAprior, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbosePlan=True, verbose=True)[1]
        o_p.writePlanToFile(oplan, filename)
    #}

    unifplan =  o_p.makeUniformExpPlan(xstart, xend, N)

    #получаем измерения с планов
    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve )


    measdataUnif = o_p.makeMeasAccToPlan_lognorm(funcf, unifplan, btrue, c,Ve )

    filename = os.path.basename(__file__).replace('.py','_results')

    #выполняем оценку
    print ("performing aprior plan")
    gknu=o_e.grandCountGN_UltraX_Qualitat (funcf, jacf,  measdata, binit, c, Ve, NSIG=100, implicit=True)
    o_q.analyseDifList(gknu, imagename='Aprior_Plan')
    o_q.printGKNUNeat(gknu)
    o_q.printQualitatNeat(measdata, gknu['b'], Ve, funcf, c, jacf)


    print ("performing uniform plan")
    gknu=o_e.grandCountGN_UltraX_Qualitat (funcf, jacf,  measdataUnif, binit, c, Ve, NSIG=100, implicit=True)
    o_q.analyseDifList(gknu, imagename='Uniform_Plan')
    o_q.printGKNUNeat(gknu)
    o_q.printQualitatNeat(measdata, gknu['b'], Ve, funcf, c, jacf)


extraction_DiodeV4_partIn_Decimal()
