

__author__ = 'reiner'

import math

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
import mpmath as mpm
from mpmath.calculus.optimization import MDNewton

import Ofiura.Ofiura_EstimationLimitedMpmath as o_elm
import Ofiura.Ofiura_planningMpmath as o_pmpm
import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Qualitat as o_q

#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit
#Опция: MPMATH USED

#Таблица переменных:

#Установка точности для mpmath. В проекте принято 50
mpm.mp.dps=40
mpm.pretty = True


#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=mpm.mpf('0.02586419') #подогнанное по pspice
foldername='cachedPlans'



def func_Kirch_DiodeV2Mod2DirectBranchMpmath(y,x,b,c):
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
    Ifwd_direct=(b[0]*(mpm.exp( (x[0]-y[0]*b[6])/(FT*b[1]))-1 ))+ (b[2]*(mpm.exp((x[0]-y[0]*b[6])/(FT*b[3]))-1 )) * (((1- (x[0]-y[0]*b[6])/b[4])**2+0.005 )**(b[5]/2))-y[0]


    return [Ifwd_direct]




def solver_Kirch_DiodeV2Mod2DirectBranchMpmath (x,b,c=None):
    """
    [Реестровая] +
    двухступенчатый - сначала np solver, затем mpm solver
    :param x:
    :param b:
    :param c:
    :return:
    """
    global numnone
    global FT


    dfdy = lambda y, x, b, c=None:  [b[2]*b[5]*b[6]*(1 - (-b[6]*y[0] + x[0])/b[4])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*(math.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)) - 1 - b[0]*b[6]*math.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + 0.005)**(b[5]/2)*math.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]
    #func = lambda y, x, b, c: [np.float(x) for x in func_Kirch_DiodeV2Mod2DirectBranchMpmath(y,x,b,c)]
    func = func_Kirch_DiodeV2Mod2DirectBranchMpmath

    solvinit=[0.5]
    try:
        solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        #solx=optimize.root(func, solvinit, args=(x,b,c),  method='lm').x
                                                                                    #solx=optimize.root(func, solvinit, args=(x,b,c), method='lm').x
                                                                                    #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
                                                                                    #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        numnone+=1
        print ("solver: ERROR first stage: "+e.__str__())
        return [None]
    funcMPM = lambda y:  func_Kirch_DiodeV2Mod2DirectBranchMpmath ([y],x,b,c)
    #dfdyMPM=lambda y: [b[2]*b[5]*b[6]*(1 - (-b[6]*[y][0] + x[0])/b[4])*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))) - 1 - b[0]*b[6]*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]
    dfdyMPM=lambda y: [b[2]*b[5]*b[6]*(1 - (-b[6]*[y][0] + x[0])/b[4])*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))) - 1 - b[0]*b[6]*mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[3]))/(FT*b[3])]
    solvinitMPM=solx[0]
    try:
        precsolx=mpm.calculus.optimization.findroot(mpm.mp, f=funcMPM, x0=solvinitMPM, solver=MDNewton, multidimensional=True, J=dfdyMPM, verify=False)
        #   precsolx=mpm.calculus.optimization.findroot(mpm.mp, f=funcMPM, x0=solvinitMPM, solver=MDNewton, multidimensional=True, verify=False)
    except BaseException as e:
        numnone+=1
        print ("solver MPM: ERROR second stage: "+e.__str__())
        return [None]
    return precsolx



def Jac_Kirch_DiodeV2Mod2DirectBranchMpmath (x,b,c,y):
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
    dfdyMPM=lambda y,x,b,c=None: mpm.matrix([[b[2]*b[5]*b[6]*(1 - (-b[6]*y[0] + x[0])/b[4])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))) - 1 - b[0]*b[6]*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]])
    dfdbMPM=lambda y,x,b,c: mpm.matrix([[mpm.exp((-b[6]*y[0] + x[0])/(FT*b[1])) - 1,
                                         -b[0]*(-b[6]*y[0] + x[0])*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                         ((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1),
                                         -b[2]*(-b[6]*y[0] + x[0])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3]**2),
                                         b[2]*b[5]*(1 - (-b[6]*y[0] + x[0])/b[4])*(-b[6]*y[0] + x[0])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]**2*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))), b[2]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)*math.log((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))/2, b[2]*b[5]*y[0]*(1 - (-b[6]*y[0] + x[0])/b[4])*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))) - b[0]*y[0]*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*y[0]*((1 - (-b[6]*y[0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*mpm.exp((-b[6]*y[0] + x[0])/(FT*b[3]))/(FT*b[3])]])
    return dfdyMPM(y,x,b,c)**-1 * dfdbMPM(y,x,b,c)


def test_Kirch_DiodeV2Mod2DirectBranchMpmath():
    btrue=[mpm.mpf('1.238e-14'), mpm.mpf('1.8'),  mpm.mpf('1.123e-14'), mpm.mpf('1.5'), mpm.mpf('1.1'), mpm.mpf('0.5'), mpm.mpf('123.')]
    b=btrue
    rng=mpm.arange('0.01','2','0.01')
    Ve=np.array([ [0.0001] ]  )
    c=None

    #снимем ВАХ
    resrng=[solver_Kirch_DiodeV2Mod2DirectBranchMpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    resrng=[o_pmpm.makeMeasOneDot_lognorm(solver_Kirch_DiodeV2Mod2DirectBranchMpmath, [x],b,c,Ve) for x in rng]
    resrng1=[solver_Kirch_DiodeV2Mod2DirectBranchMpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    plt.plot(rng , resrng, label='r=1000')
    plt.plot(rng , resrng1, label='r=3000')
    plt.legend(loc='upper left')
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()



def extraction_Kirch_DiodeV2Mod2DirectBranchMpmath():
    """
    [Реестровая]

    :return:
    """
    global FT
    global foldername

    funcf=solver_Kirch_DiodeV2Mod2DirectBranchMpmath
    jacf = Jac_Kirch_DiodeV2Mod2DirectBranchMpmath
    c=None
    Ve=np.array([ [1e-5] ]  )
    #       0           1       2     3   4    5    6
    btrue=[mpm.mpf('1.238e-14'), mpm.mpf('1.8'),  mpm.mpf('1.123e-14'), mpm.mpf('1.5'), mpm.mpf('1.1'), mpm.mpf('0.852'), mpm.mpf('123.')]
    bstart = [item-item*.2 for item in btrue]
    bend = [item+item*.22 for item in btrue]
    #binit=[mpm.mpf('1e-14'), mpm.mpf('1.1'),  mpm.mpf('1e-14'), mpm.mpf('1.'), mpm.mpf('1.'), mpm.mpf('0.5'), mpm.mpf('100.')]
    binit=mpm.matrix([['1e-14', '1.1',  '1e-14', '1.', '1.', '0.5', '100.']]).T

    xstart=[mpm.mpf('0.0001')]
    xend=[mpm.mpf('1.3')]

    N=40 #число точек в плане (для планов, кроме априорного)
    NArprior=20 #число точек в априорном плане



    #Получаем априорный план
    # import os
    # filename =foldername+'/'+os.path.basename(__file__).replace('.py','N{0}_plan'.format (NArprior))
    # print (filename)
    # try:
    #     oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
    #     print ("Read file successful")
    # except BaseException as e:
    #     oplan=o_ap.grandApriornPlanning (xstart, xend, NArprior, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbose=True)[1]
    #     o_p.writePlanToFile(oplan, filename)

    filename='cachedPlans/DiodeV2Mod2Part1N20_plan'
    try:
        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла - читаем  библиотечный файл!
        print ("Read file successful")
    except:
        exit(0)


    # newxstart=1.4
    # oplan = [item for item in oplan if item[0]<newxstart]
    #
    oplanmpm=o_pmpm.planToMpm(oplan)
    measdata = o_pmpm.makeMeasAccToPlan_lognorm(funcf, oplanmpm, btrue, c, Ve)
    gknuxlimmpm = o_elm.grandCountGN_UltraX1_Limited_wrapperMpmath(funcf,jacf,measdata,binit,bstart,bend, c, implicit=True, verbose=False, verbose_wrapper=True, isBinitGood=True )
    print (gknuxlimmpm)
    gknuxlimmpm2=o_q.convertToQualitatStandart (gknuxlimmpm, funcf, jacf,  measdata, c, Ve, name='Limited Count Aprior')
    o_q.printQualitatStandart (gknuxlimmpm2)


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
extraction_Kirch_DiodeV2Mod2DirectBranchMpmath()

