__author__ = 'reiner'
import math

import matplotlib.pyplot as plt
import numpy as np

import Ofiura_Estimation as o_e
import Ofiura_EstimationMpmath as o_empm
import Ofiura_planning as o_p
import Ofiura_ApriorPlanning as o_ap
import Ofiura_Qualitat as o_q

from scipy import optimize





numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


"""
Экстрагируем три параметра диода: N (коэфф. неидеальности), Is, R(омическое сопротивление, включенное последовательно с диодом)

"""

#Внимание - надо исключать из плана эксперимента неверно оценённые точки. Неверными признаются те, у которых результат совпадает с init
#В этом случае рабочая функция должна возвращать none, а создатель плана - вежливо вытряхивать точку из measdata


def func(y,x,b,c):
    FT=0.02586419 #подогнанное по pspice



    mm=float(b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0])

    return [mm]


resultsOfEstimation=dict() #глобальный словарь всех результатов вычислений. Из него выбирается начальное приблжение
#вектор x к вектору

def diodeResistorIMPLICITfunction (x,b,c=None):
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

    # V=x[0] #напряжение на диоде
    # Is=b[0]
    # N=b[1]
    # R=b[2]
    FT=0.02586419 #подогнанное по pspice

    dfdy=lambda y,x,b,c=None: np.array ([[ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]])
    #function=lambda y,x,b,c=None: [b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0]] #в подготовленном для взятия производной виде



    solvinit=[1]

    try:
        solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        print ('Error in findroot=',e)
        numnone+=1
        return None

    if solx-solvinit==[0]*len(solx):
         numnone+=1
         return None

    return solx

#funstr = "b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0" #в подготовленном для взятия производной виде

#print (sympy.diff(funstr,'y0').__str__())


#Производные по b
#exp((-b2*y0 + x0)/(FT*b1)) - 1
#-b0*(-b2*y0 + x0)*exp((-b2*y0 + x0)/(FT*b1))/(FT*b1**2)
#-b0*y0*exp((-b2*y0 + x0)/(FT*b1))/(FT*b1)


#y0
#-1 - b0*b2*exp((-b2*y0 + x0)/(FT*b1))/(FT*b1)


def diodeResistorIMPLICITJac (x,b,c,y):
    FT=0.02586419 #подогнанное по pspice

    dfdb=lambda y,x,b,c: np.matrix( [ [math.exp((-b[2]*y[0] + x[0])/(FT*b[1])) - 1,
                                      -b[0]*(-b[2]*y[0] + x[0])*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                      -b[0]*y[0]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]     ])


    dfdy=lambda y,x,b,c: np.matrix ([ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])])

    #возвращает структурную матрицу
    #jacf=lambda x,b,c,y: jjacf(x,b,c,y,dfdb,dfdy)


    jacf=np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


    return jacf


def testDiodeParameterExtractionIMPLICIT():
    """
    пробуем экстрагировать коэффициенты из модели диода
    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, параллельное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода
    :return:
    """


    #возвращает значение y
    funcf=diodeResistorIMPLICITfunction
    jacf = diodeResistorIMPLICITJac

    #теперь попробуем сделать эксперимент.
    c={}
    Ve=np.array([ [0.00000001] ]  )

    Ve=None

    btrue=[1.238e-14, 1.8, 100]


    bstart=np.array(btrue)-np.array(btrue)*0.2
    bend=np.array(btrue)+np.array(btrue)*0.2
    binit=[1.1e-10, 1.1,50 ]

    xstart=[0.01]
    #xend=[20,60]
    xend=[1]

    N=60
    print("performing normal research:")

    global numnone
    numnone=0
    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)

    print('unsuccessful estimations: ',numnone)



    planplot1=[x[0] for x in startplan]
    measplot1=[x['y'][0] for x in measdata]

    print(measdata)
    print(planplot1 )
    print(measplot1)


    plt.plot(planplot1, measplot1,  'ro')


    plt.ylabel('value')
    plt.xlabel('point')
    plt.grid()
    plt.show()


    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6, sign=0)
    #как мы помним, в случае неявных функций должно ставить sign=0

    print (gknu[0])

    #plotting Sk graph
    #TODO better organize: this code to estimation or somewhere
    rng=np.arange(0,len(gknu[3]))
    plt.plot(rng , gknu[3], label='Sk drop')
    plt.legend(loc='upper left')
    plt.ylabel('Sk')
    plt.xlabel('Interation')
    plt.grid()
    plt.show()






    exit(0)

    N=10
    print("performing aprior plan:")
    #oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6)[1]
    #o_p.writePlanToFile(oplan)
    oplan=o_p.readPlanFromFile() #переключение на чтение априорного плана из файла
    measdata = o_p.makeMeasAccToPlan(funcf, oplan, btrue, c,Ve )
    planplot=[x[0] for x in oplan]
    measplot=[x['y'][0] for x in measdata]



    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    measdata1 = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)








    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6, sign=1)
    print (gknu[0])
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))


        #plotting Sk graph
    #TODO better organize: this code to estimation or somewhere
    rng=np.arange(0,len(gknu[3]))
    plt.plot(rng , gknu[3], label='Sk drop')
    plt.legend(loc='upper left')
    plt.ylabel('Sk')
    plt.xlabel('Interation')




def testDiodeImplicit():

    b=[1.238e-14, 1.8, 1000]
    #0.0026638081177255196
    rng=np.arange(0.01,1,0.01)
    #снимем ВАХ
    resrng=[diodeResistorIMPLICITfunction ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
#    resrngorig=[casesDiode.diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.


    b=[1.238e-14, 1.8, 3000]


    resrng1=[diodeResistorIMPLICITfunction ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

 #   plt.plot(rng , resrngorig, label='r=0')

    plt.plot(rng , resrng, label='r=1000')
    plt.plot(rng , resrng1, label='r=3000')




    plt.legend(loc='upper left')

    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()


#тест модели
#testDiodeImplicit()

testDiodeParameterExtractionIMPLICIT()

def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)

