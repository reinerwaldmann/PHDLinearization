__author__ = 'vasilev_is'

import math
import os

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import Ofiura.Ofiura_Estimation as o_e
import Ofiura.Ofiura_ApriorPlanning as o_ap
import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Qualitat as o_q



#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit
#Таблица переменных:



#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=0.02586419 #подогнанное по pspice




def func_Kirch_DiodeV2Mod2DirectBranch(y,x,b,c):
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
    In =   float(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )) #+
    King = float(math.sqrt (b[2]/(b[2]+In))) #+
    Irec = float(b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )) #+
    Kgen = float(((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2) ) #+
    Ifwd = float(In*King+Irec*Kgen - y [0])  #+
    Ifwd_direct=(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )) * math.sqrt (b[2]/(b[2]+(b[0]*math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1)))+\
                (b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )) * (((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2))-y[0]

    return [Ifwd]


def solver_Kirch_DiodeV2Mod2DirectBranch (x,b,c=None):
    """
    [Реестровая] +
    :param x:
    :param b:
    :param c:
    :return:
    """
    global numnone

    global FT

    dfdy=lambda y,x,b,c=None: np.array ([[b[3]*b[6]*b[7]*(1 - (-b[7]*y[0] + x[0])/b[5])*
    ((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*(math.exp((-b[7]*y[0] +
    x[0])/(FT*b[4])) - 1)/(b[5]*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)) - 1 +
    b[0]**2*b[7]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*
    (math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1)*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(2*FT*b[1]*(b[0]*
    math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)) - b[0]*b[7]*
    math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*
    math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[3]*b[7]*
    ((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*math.exp((-b[7]*y[0] + x[0])/(FT*b[4]))/(FT*b[4])]])

    func = func_Kirch_DiodeV2Mod2DirectBranch
    solvinit=[1]

    try:
        solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        #print ('diodeResistorIMPLICITfunction: Error in findroot=',e)
        numnone+=1
        print ("solver: e"+e)
        return None

    if solx-solvinit==[0]*len(solx):
         numnone+=1
         print ("solver: problems in estimation")

         return None

    if solx<0:
        print ("solver: <0")
        return None

    return solx



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


    dfdb=lambda y,x,b,c: np.matrix( [[-b[0]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*(math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1)*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(2*(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)) + math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*(math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1),
                                      b[0]**2*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*(-b[7]*y[0] + x[0])*(math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1)*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(2*FT*b[1]**2*(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)) - b[0]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*(-b[7]*y[0] + x[0])*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                      b[0]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*(-b[2]/(2*(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)**2) + 1/(2*(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)))*(math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1)*(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)/b[2],
                                      ((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*(math.exp((-b[7]*y[0] + x[0])/(FT*b[4])) - 1),
                                      -b[3]*(-b[7]*y[0] + x[0])*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*math.exp((-b[7]*y[0] + x[0])/(FT*b[4]))/(FT*b[4]**2),
                                      b[3]*b[6]*(1 - (-b[7]*y[0] + x[0])/b[5])*(-b[7]*y[0] + x[0])*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*(math.exp((-b[7]*y[0] + x[0])/(FT*b[4])) - 1)/(b[5]**2*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)),
                                      b[3]*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*(math.exp((-b[7]*y[0] + x[0])/(FT*b[4])) - 1)*math.log((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)/2,
                                      b[3]*b[6]*y[0]*(1 - (-b[7]*y[0] + x[0])/b[5])*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*(math.exp((-b[7]*y[0] + x[0])/(FT*b[4])) - 1)/(b[5]*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)) + b[0]**2*y[0]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*(math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1)*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(2*FT*b[1]*(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)) - b[0]*y[0]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[3]*y[0]*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*math.exp((-b[7]*y[0] + x[0])/(FT*b[4]))/(FT*b[4])]]
   )


    dfdy=lambda y,x,b,c=None: np.array ([[b[3]*b[6]*b[7]*(1 - (-b[7]*y[0] + x[0])/b[5])*
    ((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*(math.exp((-b[7]*y[0] +
    x[0])/(FT*b[4])) - 1)/(b[5]*((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)) - 1 +
    b[0]**2*b[7]*math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*
    (math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) - 1)*math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(2*FT*b[1]*(b[0]*
    math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1)) - b[0]*b[7]*
    math.sqrt(b[2]/(b[0]*math.exp((-b[7]*y[0] + x[0])/(FT*b[1])) + b[2] - 1))*
    math.exp((-b[7]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[3]*b[7]*
    ((1 - (-b[7]*y[0] + x[0])/b[5])**2 + 0.005)**(b[6]/2)*math.exp((-b[7]*y[0] + x[0])/(FT*b[4]))/(FT*b[4])]])
    #возвращает структурную матрицу
    #jacf=lambda x,b,c,y: jjacf(x,b,c,y,dfdb,dfdy)
    jacf=np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


    return jacf


def test_Kirch_DiodeV2Mod2DirectBranch():
    #       0       1   2      3    4  5   6    7
    b=[1.238e-14, 1.8, 10E5, 1.1e-14, 2, 1, 0.5, 100]
    rng=np.arange(0.01,1.5,0.01)
    #снимем ВАХ
#    resrng=[solver_Kirch_DiodeV2Mod2DirectBranch ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    resrng=[solver_Kirch_DiodeV2Mod2DirectBranch ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

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
    funcf=solver_Kirch_DiodeV2Mod2DirectBranch
    jacf = Jac_Kirch_DiodeV2Mod2DirectBranch
    c={}
    Ve=np.array([ [0.00000000001] ]  )
    btrue=[1.238e-14, 1.8, 10E5, 1.1e-14, 2, 1, 0.5, 100]
    bstart=np.array(btrue)-np.array(btrue)*0.3
    bend=np.array(btrue)+np.array(btrue)*0.3
    binit=[1.1-14, 1.5, 9.5E5, 1e-14, 1.9, 1.1, 0.5, 100]
    xstart=[0.01]
    #xend=[20,60]
    xend=[1.1]
    N=40 #число точек в плане (для планов, кроме априорного)
    NArprior=20 #число точек в априорном плане


    #Получаем априорный план
    print("performing aprior plan:")

    #блок кеширования априорного плана в файл
    filename = os.path.basename(__file__).replace('.py','_plan')


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
    #
    # for i in range(30):
    #     measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve)
    #     gknu=o_e.grandCountGN_UltraX_ExtraStart (funcf, jacf,  measdata, bstart, bend, c, Ve,  NSIG=100, implicit=True, verbose=False, Ntries=10, name='aprior plan plus several measurements')
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
    # o_q.analyseDifList(be)
    # o_q.printGKNUNeat(be)
    # o_q.printQualitatNeat(measdata, be[0], Ve, funcf, c)
    #


    #Оценивание с использованием binit
    print('Aprior Plan Binit')
    #данные по новому формату
    binit=[1.238e-14, 1.8, 10E5, 1.1e-14, 2, 1, 0.5, 100]

    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve)
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=100, implicit=True)
    o_q.analyseDifList(gknu)
    o_q.printGKNUNeat(gknu)
    o_q.printQualitatNeat(measdata, gknu[0], Ve, funcf, c)








def solver_Kirch_DiodeV2Mod2DirectBranch_EXPLICIT(x,b,c=None):
    """
    [Реестровая]
    Уравнение Кирхгофа
    :param y:
    :param x:
    :param b:
    :param c:
    :return:
    """
    global FT
    #mm=float(b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0])
    In =   b[0]*(math.exp(x[0]/(FT*b[1]))-1 ) #+
    King = math.sqrt (b[2]/(b[2]+In)) #+
    Irec = b[3]*(math.exp((x[0])/(FT*b[4]))-1 ) #+
    Kgen = ((1- (x[0])/b[5])**2+0.005 )**(b[6]/2)  #+
    Ifwd = In*King+Irec*Kgen  #+
    #Ifwd = In  #+

    return [Ifwd]



funstr="(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )) * math.sqrt (b[2]/(b[2]+(b[0]*math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1)))+\
                (b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )) * (((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2))-y[0]"


#test_Kirch_DiodeV2Mod2DirectBranch()
# dfdb dfdy

#lst = o_d.makeDerivMatrix([funstr],list(range(8)), 'b')

#print (lst)

#Электроды

extraction_Kirch_DiodeV2Mod2DirectBranch()