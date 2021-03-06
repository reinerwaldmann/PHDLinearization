__author__ = 'reiner'
import math

import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

import Ofiura.Ofiura_Estimation as o_e
import Ofiura.Ofiura_ApriorPlanning as o_ap
import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Qualitat as o_q


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
        #print ('diodeResistorIMPLICITfunction: Error in findroot=',e)
        numnone+=1
        return None

    if solx-solvinit==[0]*len(solx):
         numnone+=1
         return None

    return solx

#funstr = "b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0" #в подготовленном для взятия производной виде

#b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0]

#print (sympy.diff(funstr,'y0').__str__())


#Производные по b
#exp((-b2*y0 + x0)/(FT*b1)) - 1
#-b0*(-b2*y0 + x0)*exp((-b2*y0 + x0)/(FT*b1))/(FT*b1**2)
#-b0*y0*exp((-b2*y0 + x0)/(FT*b1))/(FT*b1)


#y0
#-1 - b0*b2*exp((-b2*y0 + x0)/(FT*b1))/(FT*b1)


def diodeResistorIMPLICITfunction_logwrapper (x,b,c=None):
    "deprecated"
    return list(map(math.log, diodeResistorIMPLICITfunction(x,b)))






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



def plotPlanAndMeas2D(measdata):

    planplot1=[x['x'][0] for x in measdata]
    measplot1=[x['y'][0] for x in measdata]
    plt.plot(planplot1, measplot1,  'ro')
    plt.ylabel('value')
    plt.xlabel('point')
    plt.grid()
    plt.show()


def testDiodeParameterExtractionIMPLICIT(plot=True):
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
    btrue=[1.238e-14, 1.8, 100]
    #btrue=[1.5e-14, 1.75, 150]

    bstart=np.array(btrue)-np.array(btrue)*0.3
    bend=np.array(btrue)+np.array(btrue)*0.3
    binit=[1.0e-14, 1.7, 70]
    #binit=[1.238e-14, 1.8, 100]

    xstart=[0.001]
    #xend=[20,60]
    xend=[2 ]
    N=20


    print("performing aprior plan:")

#примитивная попытка автоматизировать, риальни надо кешировать в файл под хешем параметров
    try:
        oplan=o_p.readPlanFromFile() #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:
        oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbosePlan=True, verbose=True)[1]
        o_p.writePlanToFile(oplan)

    #получаем измеренные данные
#     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve )
#     #чертим эти данные
#     #o_pl.plotPlanAndMeas2D(measdata, 'Aprior Disp{0} measdata'.format(Ve))
#
#     #оценка
#     gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=True)
#
#     #вывод данных оценки - данные, квалитат, дроп
#     o_q.printGKNUNeat(gknu)
#     o_q.printQualitatNeat(measdata, gknu[0], Ve, funcf, c)
#     if plot:
#         o_pl.plotSkGraph(gknu, 'Aprior Research Sk drop')
#
#
    terminationOptDict={'VdShelfPow':-7}
#
#
#     N=30
#
#     print("\n\nperforming sequence plan with random as seed:")
#     seqplanb=o_sp.getbSeqPlanUltra(xstart, xend, N, btrue, binit, c, Ve, jacf, funcf,  dotlim=100, verbose=False, NSIG=100, implicit=True, lognorm=True, terminationOptDict=terminationOptDict) #создаём последовательный план, с выводом инфо по итерациям
#     #выводим данные последовательного планирования
#     o_q.printSeqPlanData(seqplanb)
#     #получаем данные измерения по этому последовательному плану
#     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, seqplanb[3], btrue, c,Ve)
#     #чертим эти  данные
#     if plot:
#         o_pl.plotPlanAndMeas2D(measdata, 'Hybrid Disp{0} measdata'.format(Ve))
#     print("GKNU for plan")
#     gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=True)
#     o_q.printGKNUNeat(gknu)
#     o_q.printQualitatNeat(measdataETALON, gknu[0], Ve, funcf, c)
#     if plot:
#         o_pl.plotSkGraph(gknu, 'Hybrid bei plan Sk drop')
#
#
#
#     print("\n\nperforming sequence plan with uniform as seed:")
#     unifplan=o_p.makeUniformExpPlan(xstart, xend, N)
#     seqplanb=o_sp.getbSeqPlanUltra(xstart, xend, N, btrue, binit, c, Ve, jacf, funcf, initplan=unifplan, dotlim=100, verbose=False, NSIG=100, implicit=True, lognorm=True, terminationOptDict=terminationOptDict) #создаём последовательный план, с выводом инфо по итерациям
#     #выводим данные последовательного планирования
#     o_q.printSeqPlanData(seqplanb)
#     #получаем данные измерения по этому последовательному плану
#     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, seqplanb[3], btrue, c,Ve)
#     #чертим эти  данные
#     if plot:
#         o_pl.plotPlanAndMeas2D(measdata, 'Hybrid Disp{0} measdata'.format(Ve))
#     print("GKNU for plan")
#     gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=True)
#     o_q.printGKNUNeat(gknu)
#     o_q.printQualitatNeat(measdataETALON, gknu[0], Ve, funcf, c)
#     if plot:
#         o_pl.plotSkGraph(gknu, 'Hybrid bei plan Sk drop')
#
#




    # print("\n\nperforming sequence plan with aprior as seed (hybrid):")
    # seqplanb=o_sp.getbSeqPlanUltra(xstart, xend, N, btrue, binit, c, Ve, jacf, funcf,  dotlim=100, verbose=False, NSIG=100, implicit=True, lognorm=True, terminationOptDict=terminationOptDict) #создаём последовательный план, с выводом инфо по итерациям
    # # #выводим данные последовательного планирования


    # o_q.printSeqPlanData(seqplanb)
    #получаем данные измерения по этому последовательному плану


    # #пробуем несколько измерений произвести
    # resarr=list()
    # #несколько раз измерения
    # t=PrettyTable (['Среднее логарифма правдоподобия','Сигма логарифма правдоподобия' , 'b','Среднее остатков по модулю'])
    # for i in range(20):
    #     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve)
    #     gknu=o_e.grandCountGN_UltraX_ExtraStart (funcf, jacf,  measdata, bstart, bend, c, Ve,  NSIG=100, implicit=True, verbose=False, Ntries=10, name='aprior plan plus several measurements')
    #     if (gknu):
    #         resarr.append(gknu)
    # if resarr:
    #     for gknu in resarr:
    #         if (gknu):
    #             t.add_row([gknu['AvLogTruth'],gknu['SigmaLT'], gknu['b'], gknu['AvDif'] ])
    # be=o_e.selectBestEstim (resarr)
    # t.add_row(['*', '*', '*', '*' ])
    # t.add_row([be['AvLogTruth'], be['SigmaLT'], be['b'], be['AvDif'] ])
    # print(t)
    # o_q.analyseDifList(be)




    #Априорный план просто


    startplan =  o_p.makeUniformExpPlan(xstart, xend, 150)

    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c,Ve)
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=True)
    o_q.printGKNUNeat(gknu)
    o_q.printQualitatNeat(measdata, gknu[0], Ve, funcf, c, jacf)






    # measdata = o_p.makeMeasAccToPlan_lognorm(funcf, seqplanb[3], btrue, c,Ve)
    #чертим эти  данные
    #if plot:
    #    o_pl.plotPlanAndMeas2D(measdata, 'Hybrid Disp{0} measdata'.format(Ve))
    #здесь мы сравниваем, сработает ли метод уточнения результатов - выбор начального значения из пула
    # print('grandCountGN_UltraX_ExtraStart')
    #
    # #данные по новому формату
    #
    #
    #
    # print('normal')
    # print("GKNU bei plan")
    # gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=True)
    # o_q.printGKNUNeat(gknu)
    # o_q.printQualitatNeat(measdata, gknu[0], Ve, funcf, c)
    # if plot:
    #     o_pl.plotSkGraph(gknu, 'Hybrid bei plan Sk drop')




    # while (True):
    #     try:
    #         N=30
    #         print("performing aprior plan:")
    #         oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbosePlan=True)[1]
    #         o_p.writePlanToFile(oplan)
    #         #oplan=o_p.readPlanFromFile() #переключение на чтение априорного плана из файла
    #         measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve )
    #         #o_pl.plotPlanAndMeas2D(measdata, 'Aprior Disp{0} measdata'.format(Ve))
    #         gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6, sign=1)
    #         print (gknu[0])
    #         print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))
    #         o_pl.plotSkGraph(gknu, 'Aprior Research Sk drop')
    #         return
    #     except BaseException as e:
    #         pass




#начерталово графика плана измерения и графика без дисперсии, по графикам видно, что дисперсия хороша
    # planplot1=[x['x'][0] for x in measdata]
    # measplot1=[x['y'][0] for x in measdata]
    # plt.plot(planplot1, measplot1,  'ro')
    # measdata_pure = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c, None)
    # measplot2=[x['y'][0] for x in measdata_pure]
    # plt.plot(planplot1, measplot2)
    # plt.ylabel('value')
    # plt.xlabel('point')
    # plt.grid()
    # plt.show()
    #
    # exit(0)


def testDiodeImplicit():

    b=[1.238e-14, 1.3, 0    ]
    #0.0026638081177255196
    rng=np.arange(0.01,2,0.01)
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
testDiodeImplicit()

#шобы при обвале по переполнению автоматом заново запускался
#clsbeep.beep()


#testDiodeParameterExtractionIMPLICIT(plot=False)






def closest_node(node, nodes):
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    return np.argmin(dist_2)



 # N=60
    # print("performing normal research:")
    #
    # global numnone
    # numnone=0
    # #создание стартового плана
    #startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    #measdata = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c, Ve)
    #o_pl.plotPlanAndMeas2D(measdata, 'Normal Uniform Disp{0} measdata'.format(Ve))
    # print('unsuccessful estimations: ',numnone)
    # gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6, sign=0)
    # #как мы помним, в случае неявных функций должно ставить sign=0
    # print (gknu[0])
    # o_pl.plotSkGraph(gknu, 'Normal Research Sk drop')
    # print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))

    #N=30
    # print("performing aprior plan:")
    # # oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbosePlan=True)[1]
    # o_p.writePlanToFile(oplan)