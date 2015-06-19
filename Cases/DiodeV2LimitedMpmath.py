from mpmath.functions.functions import re

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
#mpm.mp.dps=40 #данные параметры уже выставлены в офиуре
#mpm.pretty = True


#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=mpm.mpf('0.02586419') #подогнанное по pspice
foldername='cachedPlans'

XO=mpm.mpf(1)
XT=mpm.mpf(2)
XF=mpm.mpf('0.005')



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


    v=x[0]
    i=y[0]

    iss=IS=b[0]
    n=N=b[1]
    ikf=IKF=b[2]
    isr=ISR=b[3]
    nr=NR=b[4]
    vj=VJ=b[5]
    m=M=b[6]
    rs=RS=b[7]

    In=IS * (math.exp((v-i*RS)/(FT*N))-1  )
    Kinj=math.sqrt( IKF/ (IKF+In) )
    Irec=ISR*(math.exp((v-i*RS)/(FT*NR))-1)
    Kgen=((1-(v-i*RS)/VJ)**2 +.005)**(M/2)
    Ifwd = In * Kinj + Irec*Kgen
    zero_sum = Ifwd -i

    First_half = IS * (math.exp((v-i*RS)/(FT*N))-1  ) *  math.sqrt( IKF/ (IKF+IS * (math.exp((v-i*RS)/(FT*N))-1  )) )
    Second_half = ISR*(math.exp((v-i*RS)/(FT*NR))-1) * ((1-(v-i*RS)/VJ)**2 +.005)**(M/2)

    zero_sum = First_half + Second_half -i

    return [zero_sum]

def dfdy (y,x,b,c):

    global FT,XO,XT,XF

    ft=FT

    v=x[0]
    i=y[0]

    iss=IS=b[0]
    n=N=b[1]
    ikf=IKF=b[2]
    isr=ISR=b[3]
    nr=NR=b[4]
    vj=VJ=b[5]
    m=M=b[6]
    rs=RS=b[7]

    #fh = iss**mpm.mpf(2)*rs*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - mpm.mpf(1))))*(mpm.exp((-i*rs + v)/(ft*n)) - mpm.mpf(1))*mpm.exp((-i*rs + v)/(ft*n))/(mpm.mpf(2)*ft*n*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - mpm.mpf(1)))) - iss*rs*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - mpm.mpf(1))))*mpm.exp((-i*rs + v)/(ft*n))/(ft*n)
    #sh = isr*m*rs*(mpm.mpf(1) - (-i*rs + v)/vj)*((mpm.mpf(1) - (-i*rs + v)/vj)**mpm.mpf(2) + mpm.mpf('0.005'))**(m/mpm.mpf(2))*(mpm.exp((-i*rs + v)/(ft*nr)) - mpm.mpf(1))/(vj*((mpm.mpf(1) - (-i*rs + v)/vj)**mpm.mpf(2) + mpm.mpf('0.005'))) - isr*rs*((mpm.mpf(1) - (-i*rs + v)/vj)**mpm.mpf(2) + mpm.mpf('0.005'))**(m/mpm.mpf(2))*mpm.exp((-i*rs + v)/(ft*nr))/(ft*nr)

    fh = iss**XT*rs*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(mpm.exp((-i*rs + v)/(ft*n)) - XO)*mpm.exp((-i*rs + v)/(ft*n))/(XT*ft*n*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))) - iss*rs*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*mpm.exp((-i*rs + v)/(ft*n))/(ft*n)
    sh = isr*m*rs*(XO - (-i*rs + v)/vj)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(mpm.exp((-i*rs + v)/(ft*nr)) - XO)/(vj*((XO - (-i*rs + v)/vj)**XT + XF)) - isr*rs*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*mpm.exp((-i*rs + v)/(ft*nr))/(ft*nr)



    return fh+sh

def dfdb (y,x,b,c):

    global FT,XO,XT,XF

    ft=FT

    v=x[0]
    i=y[0]

    iss=IS=b[0]
    n=N=b[1]
    ikf=IKF=b[2]
    isr=ISR=b[3]
    nr=NR=b[4]
    vj=VJ=b[5]
    m=M=b[6]
    rs=RS=b[7]

    #ещё раз: ВСЕ константы должны быть в перспективе либо глобальными переменными, как FT, либо составляющими вектора c, но он кривой, потому лучше уж глобальными.

    return mpm.matrix ([[iss*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(-mpm.exp((-i*rs + v)/(ft*n)) + XO)*(mpm.exp((-i*rs + v)/(ft*n)) - XO)/(XT*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))) + mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(mpm.exp((-i*rs + v)/(ft*n)) - XO),
                         iss**XT*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(-i*rs + v)*(mpm.exp((-i*rs + v)/(ft*n)) - XO)*mpm.exp((-i*rs + v)/(ft*n))/(XT*ft*n**XT*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))) - iss*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(-i*rs + v)*mpm.exp((-i*rs + v)/(ft*n))/(ft*n**XT),
                        iss*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))*(-ikf/(XT*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))**XT) + XO/(XT*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))))*(mpm.exp((-i*rs + v)/(ft*n)) - XO)/ikf,
                         ((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(mpm.exp((-i*rs + v)/(ft*nr)) - XO),
                         -isr*(-i*rs + v)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*mpm.exp((-i*rs + v)/(ft*nr))/(ft*nr**XT),
                        isr*m*(XO - (-i*rs + v)/vj)*(-i*rs + v)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(mpm.exp((-i*rs + v)/(ft*nr)) - XO)/(vj**XT*((XO - (-i*rs + v)/vj)**XT + XF)),
                         isr*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(mpm.exp((-i*rs + v)/(ft*nr)) - XO)*mpm.log((XO - (-i*rs + v)/vj)**XT + XF)/XT,
                         i*iss**XT*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*(mpm.exp((-i*rs + v)/(ft*n)) - XO)*mpm.exp((-i*rs + v)/(ft*n))/(XT*ft*n*(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO))) - i*iss*mpm.sqrt(ikf/(ikf + iss*(mpm.exp((-i*rs + v)/(ft*n)) - XO)))*mpm.exp((-i*rs + v)/(ft*n))/(ft*n)+ \
                         i*isr*m*(XO - (-i*rs + v)/vj)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(mpm.exp((-i*rs + v)/(ft*nr)) - XO)/(vj*((XO - (-i*rs + v)/vj)**XT + XF)) - i*isr*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*mpm.exp((-i*rs + v)/(ft*nr))/(ft*nr)
                        ]])

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
    return dfdy(y,x,b,c)**mpm.mpf(-1) * dfdb(y,x,b,c)

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


    func = lambda y, x, b, c:  [float(func_Kirch_DiodeV2Mod2DirectBranchMpmath (y,x,b,c)[0])]


    dfdynumpy=lambda y,x,b,c=None: np.array ([[  float( dfdy(y,x,b,c) )     ]])

    solvinit=[0.001]
    try:
        solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdynumpy, method='lm').x
        #solx=optimize.root(func, solvinit, args=(x,b,c),  method='lm').x
                                                                                    #solx=optimize.root(func, solvinit, args=(x,b,c), method='lm').x
                                                                                    #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
                                                                                    #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        numnone+=1
        print ("solver: ERROR first stage: "+e.__str__())
        return [None]
    funcMPM = lambda y:  func_Kirch_DiodeV2Mod2DirectBranchMpmath ([y],x,b,c)

    dfdyMPM=lambda y: [b[2]*b[5]*b[6]*(1 - (-b[6]*[y][0] + x[0])/b[4])*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*(mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[3])) - 1)/(b[4]*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))) - 1 - b[0]*b[6]*mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[1]))/(FT*b[1]) - b[2]*b[6]*((1 - (-b[6]*[y][0] + x[0])/b[4])**2 + mpm.mpf('0.005'))**(b[5]/2)*mpm.exp((-b[6]*[y][0] + x[0])/(FT*b[3]))/(FT*b[3])]
    #dfdyMPM=dfdy

    solvinitMPM=mpm.mpf(solx[0].__str__())
    try:
        precsolx=mpm.calculus.optimization.findroot(mpm.mp, f=funcMPM, x0=solvinitMPM, solver=MDNewton, multidimensional=True, J=dfdyMPM, verify=False)

    except BaseException as e:
        numnone+=1
        print ("solver MPM: ERROR second stage: "+e.__str__())

        return [None]
    return precsolx














def test_Kirch_DiodeV2Mod2DirectBranchMpmath():
    btrue=[mpm.mpf('341.4e-6'), mpm.mpf('2.664'), mpm.mpf('37.08e-3'), mpm.mpf('17.26e-27'), mpm.mpf('5.662'), mpm.mpf('4.282'), mpm.mpf('0.5751'), mpm.mpf('3.65e-3')]

    btrueNormalDiode=[mpm.mpf('6.233e-8'), mpm.mpf('2.1'), mpm.mpf('1000'), mpm.mpf('17.26e-100'), mpm.mpf('2'), mpm.mpf('1'), mpm.mpf('0.5'), mpm.mpf('3.665e-2')]

    #btrueNormalDiode=[mpm.mpf('7.69e-8'), mpm.mpf('1.45') , mpm.mpf('10'), mpm.mpf(0),mpm.mpf(2),mpm.mpf(1),mpm.mpf('0.5'), mpm.mpf('.0422')] #номинальные значения диода D1N4001 с сайта, вроде официальной модели производителя

    b=btrue
    rng=mpm.arange('0.001','1.5','0.01')
    Ve=np.array([ [1.9e-5] ]  )
    c=None

    #снимем ВАХ
    #resrng=[solver_Kirch_DiodeV2Mod2DirectBranchMpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

    resrng=list()
    nrng=list()

    resrng1,nrng1=list(), list()

    for x in rng:
        slv=solver_Kirch_DiodeV2Mod2DirectBranchMpmath ([x],b)

        slv1=solver_Kirch_DiodeV2Mod2DirectBranchMpmath ([x],btrueNormalDiode)

        if slv[0]:
            if (slv[0]>=1):
                pass
                #print (x)
                #break
            resrng.append(slv[0])
            nrng.append(x)

        else:
            print (slv, x,b)


        if slv1[0]:
            if (slv1[0]>=1):
                pass
                #print (x)
                #break
            resrng1.append(slv1[0])
            nrng1.append(x)


        else:
            print (slv1, x,b)

    #resrng=[o_pmpm.makeMeasOneDot_lognorm(solver_Kirch_DiodeV2Mod2DirectBranchMpmath, [x],b,c,Ve) for x in rng]


    #resrng1=[solver_Kirch_DiodeV2Mod2DirectBranchMpmath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.


    plt.plot(nrng1 , resrng1, label='normal')
    plt.plot(nrng , resrng, label='Schottky')
    plt.legend(loc='upper left')

    plt.axis([0,1,0,1])

    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()



def extraction_Kirch_DiodeV2Mod2DirectBranchMpmath():
    """
    [Реестровая]

    :return:
    """
    global FT,XO,XT,XF
    global foldername

    funcf=solver_Kirch_DiodeV2Mod2DirectBranchMpmath
    jacf = Jac_Kirch_DiodeV2Mod2DirectBranchMpmath
    c=None
    Ve=np.array([ [1e-5] ]  )
    #       0           1       2     3   4    5    6
    btrue=[mpm.mpf('341.4e-6'), mpm.mpf('2.664'), mpm.mpf('37.08e-3'), mpm.mpf('17.26e-27'), mpm.mpf('5.662'), mpm.mpf('4.282'), mpm.mpf('0.5751'), mpm.mpf('3.65e-3')]


    #коридор в 20 процентов в обе стороны
    diap = mpm.mpf('.2')
    bstart = [item-item*diap for item in btrue]
    bend = [item+item*diap for item in btrue]


    binit=btrue


    xstart=[mpm.mpf('0.001')]
    xend=[mpm.mpf('1.5')]

    N=20

    print("performing aprior plan:")

#примитивная попытка автоматизировать, риальни надо кешировать в файл под хешем параметров

    import os
    filename =foldername+'/'+'_'.format(N)+os.path.basename(__file__).replace('.py','_plan')

    try:
        oplan=o_p.readPlanFromFile(filename) #переключение на чтение априорного плана из файла
        print ("Read file successful")
    except BaseException as e:
        oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbose=True)[1]
        o_p.writePlanToFile(oplan, filename)









exit(0)

First_half = 'ISS * (math.exp((v-i*RS)/(FT*N))-1  ) *  math.sqrt( IKF/ (IKF+ISS * (math.exp((v-i*RS)/(FT*N))-1  )) )'
Second_half = 'ISR*(math.exp((v-i*RS)/(FT*NR))-1) * ((1-(v-i*RS)/VJ)**2 +.005)**(M/2)'
import sympy as smp

ret_normal = lambda  dfstr: dfstr.replace('exp', 'math.exp').replace('sqrt', 'math.sqrt').replace ('log', 'math.log')
dff = lambda str1,var:  str(smp.diff ( str1.replace('math.','').lower(), var)).replace('exp', 'math.exp').replace('sqrt', 'math.sqrt').replace ('log', 'math.log')
dffmpm = lambda str1,var:  str(smp.diff ( str1.replace('math.','').lower(), var)).replace('exp', 'math.exp').replace('sqrt', 'math.sqrt').replace ('log', 'math.log').replace('math','mpm')
varlist = ['iss', 'n','ikf','isr','nr','vj','m','rs']

for var in varlist:
    print (dffmpm(First_half,var))
    print (dffmpm(Second_half,var))
    print ()















