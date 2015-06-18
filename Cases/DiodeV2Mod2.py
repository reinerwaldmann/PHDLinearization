__author__ = 'vasilev_is'

import math
import os

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import Ofiura.Ofiura_EstimationLimited as o_el
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
    In =   float(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )) #+
    King = float(math.sqrt (b[2]/(b[2]+In))) #+
    Irec = float(b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )) #+
    Kgen = float(((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2) ) #+
    Ifwd = float(In*King+Irec*Kgen - y [0])  #+



    #модель задана неявно
    zero_sum=(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )) *\
                math.sqrt (b[2]/(b[2]+(b[0]*math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1)))+\
                (b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )) * \
                (((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2)) -y[0]

    #Is


    return [zero_sum]


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
        print ("solver: "+e.__str__())
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
    #b=[341.4e-6, 2.664, 37.08e-3, 17.26e-27, 5.662, 4.282, 0.5751, 3.65e-3]

    b=[341.4e-6, 2.664, 37.08e-3, 17.26e-27, 5.662, 4.282, 0.5751, 3.65e-3]

    rng=np.arange(0.6, .77, 0.0001)
    #снимем ВАХ
#    resrng=[solver_Kirch_DiodeV2Mod2DirectBranch ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.



    resrng = list()
    nrng = list()

    for x in rng:
        slv=solver_Kirch_DiodeV2Mod2DirectBranch ([x],b)

        if slv:
            if (slv[0]>=1):
                print (x)
                break

            resrng.append(slv[0])
            nrng.append(x)
        else:

            print (slv, x,b)

    #resrng=[solver_Kirch_DiodeV2Mod2DirectBranch ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

#    resrngorig=[casesDiode.diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    plt.plot(nrng , resrng, label='schottky 10bq100')
    plt.legend(loc='upper left')
    plt.axis([0.0,1.0,0,1])
    plt.grid()
    plt.show()


def plotPlanAndMeas2D(measdata):

    planplot1=[x['x'][0] for x in measdata]
    measplot1=[x['y'][0] for x in measdata]
    plt.plot(planplot1, measplot1,  'r')
    plt.ylabel('value')
    plt.xlabel('point')
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
    Ve=np.array([[1.9e-5]])

    btrue=[341.4e-6, 2.664, 37.08e-3, 17.26e-27, 5.662, 4.282, 0.5751, 3.65e-3]

    bstart=np.array(btrue)-np.array(btrue)*0.05
    bend=np.array(btrue)+np.array(btrue)*0.05

    binit=[341.4e-6, 2.664, 37.08e-3, 17.26e-27, 5.662, 4.282, 0.5751, 3.65e-3] #binit=btrue так как взято из спайс-модели



    xstart=[0.7]
    #xend=[20,60]
    xend=[.77]

    N=20 #число точек в априорном плане


    #Получаем априорный план
    print("performing aprior plan:")

    #блок кеширования априорного плана в файл
    import os
    filename =foldername+'/'+'RD_10BQ100_N{0}_Dev-62012P_'.format(N)+os.path.basename(__file__).replace('.py','_plan')

    try:

        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:
        oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbose=True)[1]
        o_p.writePlanToFile(oplan, filename)

    # #Задание опций для последовательного плана
    # terminationOptDict={'VdShelfPow':-7}

    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, oplan, btrue, c,Ve )
    plotPlanAndMeas2D(measdata)

    gknuxlim = o_el.grandCountGN_UltraX1_Limited_wrapper(funcf,jacf,measdata,binit,bstart,bend, c, implicit=True, verbose=False, verbose_wrapper=False )
    gknuxlim2=o_q.convertToQualitatStandart (gknuxlim, funcf, jacf,  measdata, c, Ve, name='Limited Count Aprior')
    o_q.printQualitatStandart (gknuxlim2)





# funstr="(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )) * math.sqrt (b[2]/(b[2]+(b[0]*math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1)))+\
#                 (b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )) * (((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2))-y[0]"


#test_Kirch_DiodeV2Mod2DirectBranch()
# dfdb dfdy

#lst = o_d.makeDerivMatrix([funstr],list(range(8)), 'b')

#print (lst)

#Электроды

#test_Kirch_DiodeV2Mod2DirectBranch()
extraction_Kirch_DiodeV2Mod2DirectBranch()