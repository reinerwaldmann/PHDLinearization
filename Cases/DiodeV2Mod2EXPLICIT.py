__author__ = 'vasilev_is'
__author__ = 'vasilev_is'

import math
import os

import numpy as np
import matplotlib.pyplot as plt

import Ofiura.Ofiura_Estimation as o_e
import Ofiura.Ofiura_ApriorPlanning as o_ap
import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Qualitat as o_q






#Part1: прямая ветвь ВАХ диода
#Режим явной функции, резистора нет

#Стандарт: explicit
#Таблица переменных:



#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=0.02586419 #подогнанное по pspice




def func_Kirch_DiodeV2Mod2DirectBranch(x,b,c=None):
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
    #mm=float(b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0])
    # In =   float(b[0]*(math.exp( (x[0]-y[0]*b[6])/(FT*b[1]))-1 )) #+
    # King = 1
    # Irec = float(b[2]*(math.exp((x[0]-y[0]*b[6])/(FT*b[3]))-1 )) #+
    # Kgen = float(((1- (x[0]-y[0]*b[6])/b[4])**2+0.005 )**(b[5]/2) ) #+
    # Ifwd = float(In*King+Irec*Kgen - y [0])  #+

    In=b[0]*(math.exp( x[0]/(FT*b[1]))-1 )
    King=1
    Irec=b[2]*(math.exp(x[0]/(FT*b[3]))-1 )
    Kgen = ((1- (x[0])/b[4])**2+0.005 )**(b[5]/2)


    print (In, King, Irec, Kgen, x,b)

    try:
        Ifwd_direct=b[0]*(math.exp( x[0]/(FT*b[1]))-1 )+ b[2]*(math.exp(x[0]/(FT*b[3]))-1 ) * (((1- (x[0])/b[4])**2+0.005 )**(b[5]/2))
    except BaseException as e:
        return None

    return [Ifwd_direct]




def Jac_Kirch_DiodeV2Mod2DirectBranch (x,b,c,y):
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

    dfdb=lambda y,x,b,c: np.matrix([[
                                        math.exp(x[0]/(FT*b[1])) - 1,
                                        -b[0]*x[0]*math.exp(x[0]/(FT*b[1]))/(FT*b[1]**2),
                                        ((1 - x[0]/b[4])**2 + 0.005)**(b[5]/2)*(math.exp(x[0]/(FT*b[3])) - 1),
                                        -b[2]*x[0]*((1 - x[0]/b[4])**2 + 0.005)**(b[5]/2)*math.exp(x[0]/(FT*b[3]))/(FT*b[3]**2),
                                        b[2]*b[5]*x[0]*(1 - x[0]/b[4])*((1 - x[0]/b[4])**2 + 0.005)**(b[5]/2)*(math.exp(x[0]/(FT*b[3])) - 1)/(b[4]**2*((1 - x[0]/b[4])**2 + 0.005)),
                                        b[2]*((1 - x[0]/b[4])**2 + 0.005)**(b[5]/2)*(math.exp(x[0]/(FT*b[3])) - 1)*math.log((1 - x[0]/b[4])**2 + 0.005)/2
                                    ]])

    #print (dfdb(y,x,b,c))

    func_Kirch_DiodeV2Mod2DirectBranch(x,b,c)


    try:
        return dfdb(y,x,b,c)
    except BaseException:
        return None


def test_Kirch_DiodeV2Mod2DirectBranch():
    #       0       1   2      3    4  5   6    7
    b=[1.238e-14, 1.8, 1.1e-20, 2.0, 1.0, 0.5]
    rng=np.arange(0.01,1.5,0.01)
    #снимем ВАХ
#    resrng=[solver_Kirch_DiodeV2Mod2DirectBranch ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    resrng=[func_Kirch_DiodeV2Mod2DirectBranch ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

#    resrngorig=[casesDiode.diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    plt.plot(rng , resrng, label='r=1000')
    plt.legend(loc='upper left')
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()



def extraction_Kirch_DiodeV2Mod2DirectBranch():
    """
    [Реестровая]

    :return:
    """
    funcf=func_Kirch_DiodeV2Mod2DirectBranch
    jacf = Jac_Kirch_DiodeV2Mod2DirectBranch
    c={}
    Ve=np.array([ [0.0000000000001] ]  )
    btrue=[1.238e-14, 1.8, 1.5e-15, 1.9, 1.2, 0.5]

    bstart=np.array(btrue)-np.array(btrue)*0.1
    bend=np.array(btrue)+np.array(btrue)*0.1
    binit=[1.238e-14, 1.7,  1.1e-20, 1.9, 1.1, 0.5]


    xstart=[0.01]
    #xend=[20,60]
    xend=[1.5]
    N=200 #число точек в плане (для планов, кроме априорного)
    NArprior=30 #число точек в априорном плане




    #Получаем априорный план
    print("performing aprior plan:")
    #блок кеширования априорного плана в файл
    filename = os.path.basename(__file__).replace('.py','_plan1')
    try:

        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:
        oplan=o_ap.grandApriornPlanning (xstart, xend, NArprior, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbosePlan=True, verbose=True)[1]
        o_p.writePlanToFile(oplan, filename)

    #Задание опций для последовательного плана
    terminationOptDict={'VdShelfPow':-7}



    #Оценивание параметров с использованием  ряда начальных значений
    # resarr=list() #Список результатов
    # t=PrettyTable (['Среднее логарифма правдоподобия','Сигма логарифма правдоподобия' , 'b','Среднее остатков по модулю'])
    # for i in range(30):
    #     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve)
    #     gknu=o_e.grandCountGN_UltraX_ExtraStart (funcf, jacf,  measdata, bstart, bend, c, Ve,  NSIG=100, implicit=False, verbose=False, Ntries=10, name='aprior plan plus several measurements')
    #     if (gknu):
    #         resarr.append(gknu)
    # if resarr:
    #     for gknu in resarr:
    #         if (gknu):
    #             t.add_row([gknu['AvLogTruth'],gknu['SigmaLT'], gknu['b'], gknu['AvDif'] ])
    # be=o_e.selectBestEstim (resarr)
    #
    # t.add_row(['*', '*', '*', '*' ])
    # t.add_row([be['AvLogTruth'], be['SigmaLT'], be['b'], be['AvDif'] ])
    # print(t)
    # #o_q.analyseDifList(be)
    # o_q.printGKNUNeat(be)
    # o_q.printQualitatNeat(measdata, be[0], Ve, funcf, c)





    #Оценивание с использованием binit

#    По априорному плану
#     print('Aprior Plan Binit')
#     #данные по новому формату
#
#     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve)
#     gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=False)
#     #o_q.analyseDifList(gknu)
#     o_q.printGKNUNeat(gknu)
#     o_q.printQualitatNeat(measdata, gknu[0], Ve, funcf, c)


    #
    # #По нормальному плану

    print("performing normal research:")
    startplan =  o_p.makeUniformExpPlan(xstart, xend, 150)
    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c,Ve)
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=False)
    if (gknu):
        #o_q.analyseDifList(gknu)
        o_q.printGKNUNeat(gknu)
        o_q.printQualitatNeat(measdata, gknu[0], Ve, funcf, c, jacf)





#funstr="(b[0]*(math.exp( x[0]/(FT*b[1]))-1 ))+ (b[2]*(math.exp(x[0]/(FT*b[3]))-1 )) * (((1- x[0]/b[4])**2+0.005 )**(b[5]/2))"


#test_Kirch_DiodeV2Mod2DirectBranch()

# dfdb dfdy
#

#
#print ("dfdb")
#print (o_d.makeDerivMatrix([funstr],list(range(6)), 'b'))


#exit(0)

extraction_Kirch_DiodeV2Mod2DirectBranch()