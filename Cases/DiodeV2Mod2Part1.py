__author__ = 'vasilev_is'

import math
import os

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import Ofiura.Ofiura_EstimationLimited as o_el
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
foldername='cachedPlans'




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
    # In =   float(b[0]*(math.exp( (x[0]-y[0]*b[6])/(FT*b[1]))-1 )) #+
    # King = 1
    # Irec = float(b[2]*(math.exp((x[0]-y[0]*b[6])/(FT*b[3]))-1 )) #+
    # Kgen = float(((1- (x[0]-y[0]*b[6])/b[4])**2+0.005 )**(b[5]/2) ) #+
    # Ifwd = float(In*King+Irec*Kgen - y [0])  #+
    Ifwd_direct=(b[0]*(math.exp( (x[0]-y[0]*b[6])/(FT*b[1]))-1 ))+ (b[2]*(math.exp((x[0]-y[0]*b[6])/(FT*b[3]))-1 )) * (((1- (x[0]-y[0]*b[6])/b[4])**2+0.005 )**(b[5]/2))-y[0]

    return [Ifwd_direct]


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

    dfdy=lambda y,x,b,c=None: np.array ([[b[2]*b[5]*b[6]*(1 - (-b[6]*y[0] + x[0])/b[4])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)) - 1 - b[0]*b[6]*math.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*math.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]]
)

    func = func_Kirch_DiodeV2Mod2DirectBranch
    solvinit=[1]

    try:
        solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        #solx=optimize.root(func, solvinit, args=(x,b,c), method='lm').x
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
        #print ("solver: <0")
        pass
        #return None

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


    dfdy=lambda y,x,b,c=None: np.array ([[b[2]*b[5]*b[6]*(1 - (-b[6]*y[0] + x[0])/b[4])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)) - 1 - b[0]*b[6]*math.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*math.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]]
    )


    dfdb=lambda y,x,b,c: np.matrix([[math.exp((-b[6]*y[0] + x[0])/(FT*b[1])) - 1, -b[0]*(-b[6]*y[0] + x[0])*math.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2), ((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1), -b[2]*(-b[6]*y[0] + x[0])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*math.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3]**2), b[2]*b[5]*(1 - (-b[6]*y[0] + x[0])/b[4])*(-b[6]*y[0] + x[0])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]**2*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)), b[2]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)*math.log((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)/2, b[2]*b[5]*y[0]*(1 - (-b[6]*y[0] + x[0])/b[4])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)) - b[0]*y[0]*math.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*y[0]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*math.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]]
)





    #возвращает структурную матрицу
    #jacf=lambda x,b,c,y: jjacf(x,b,c,y,dfdb,dfdy)
    jacf=np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


    return jacf


def test_Kirch_DiodeV2Mod2DirectBranch():
    #       0       1   2      3    4  5   6    7
    btrue=[1.238e-14, 1.8,  1.123e-14, 1.5, 1., 0.5, 100.]

    bstart=np.array(btrue)-np.array(btrue)*0.1
    bend=np.array(btrue)+np.array(btrue)*0.102

    b=btrue

    rng=np.arange(0.01,1.5,0.001)
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
    global FT
    global foldername

    funcf=solver_Kirch_DiodeV2Mod2DirectBranch
    jacf = Jac_Kirch_DiodeV2Mod2DirectBranch
    c={}
    Ve=np.array([ [0.00000000001] ]  )
    #       0           1       2     3   4    5    6
    btrue=[1.238e-14, 1.8,  1.123e-14, 1.5, 1., 0.5, 100.]


    bstart=np.array(btrue)-np.array(btrue)*0.1
    bend=np.array(btrue)+np.array(btrue)*0.102
    binit=[1.1-14, 1.5,  1.1e-14, 1.9, 1.1, 0.8, 98.]
    xstart=[0.001]
    xend=[1.5]
    N=40 #число точек в плане (для планов, кроме априорного)
    NArprior=20 #число точек в априорном плане


    #Получаем априорный план
    import os
    filename =foldername+'/'+os.path.basename(__file__).replace('.py','_plan')
    try:
        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:
        oplan=o_ap.grandApriornPlanning (xstart, xend, NArprior, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbose=True)[1]
        o_p.writePlanToFile(oplan, filename)

    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve )


    gknuxlim = o_el.grandCountGN_UltraX1_Limited_wrapper(funcf,jacf,measdata,binit,bstart,bend, c, implicit=True, verbose=False, verbose_wrapper=False )
    gknux = o_e.grandCountGN_UltraX1(funcf, jacf, measdata, binit, c, implicit=True)

    gknuxlim2=o_q.convertToQualitatStandart (gknuxlim, funcf, jacf,  measdata, c, Ve, name='Limited Count Aprior')
    gknux2=o_q.convertToQualitatStandart (gknux, funcf, jacf,  measdata, c, Ve, name='Normal Count Aprior')

    o_q.printQualitatStandart (gknuxlim2)
    o_q.printQualitatStandart (gknux2)





funstr="(b[0]*(math.exp( (x[0]-y[0]*b[6])/(FT*b[1]))-1 ))+ (b[2]*(math.exp((x[0]-y[0]*b[6])/(FT*b[3]))-1 )) * (((1- (x[0]-y[0]*b[6])/b[4])**2+0.005 )**(b[5]/2))-y[0]"


#test_Kirch_DiodeV2Mod2DirectBranch()
# dfdb dfdy
#
# print ("dfdy")
# print (o_d.makeDerivMatrix([funstr],[0], 'y'))
#
# print ("dfdb")
# print (o_d.makeDerivMatrix([funstr],list(range(8)), 'b'))
# exit(0)
#

#test_Kirch_DiodeV2Mod2DirectBranch()
extraction_Kirch_DiodeV2Mod2DirectBranch()