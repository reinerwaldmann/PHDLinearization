
__author__ = 'reiner'

import math

import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

import Ofiura.Ofiura_ApriorPlanning as o_ap

#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit
#Опция: mathATH USED

#Таблица переменных:

#Установка точности для mathath. В проекте принято 50
#math.mp.dps=40 #данные параметры уже выставлены в офиуре
#math.pretty = True


#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=0.02586419 #подогнанное по pspice
foldername='cachedPlans'

XO=1
XT=2
XF=0.005

 
def func_Kirch_DiodeV2Mod2DirectBranchmathath(y,x,b,c):
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

    print (First_half, Second_half, i)
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

    #fh = iss**math.mpf(2)*rs*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - math.mpf(1))))*(math.exp((-i*rs + v)/(ft*n)) - math.mpf(1))*math.exp((-i*rs + v)/(ft*n))/(math.mpf(2)*ft*n*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - math.mpf(1)))) - iss*rs*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - math.mpf(1))))*math.exp((-i*rs + v)/(ft*n))/(ft*n)
    #sh = isr*m*rs*(math.mpf(1) - (-i*rs + v)/vj)*((math.mpf(1) - (-i*rs + v)/vj)**math.mpf(2) + math.mpf('0.005'))**(m/math.mpf(2))*(math.exp((-i*rs + v)/(ft*nr)) - math.mpf(1))/(vj*((math.mpf(1) - (-i*rs + v)/vj)**math.mpf(2) + math.mpf('0.005'))) - isr*rs*((math.mpf(1) - (-i*rs + v)/vj)**math.mpf(2) + math.mpf('0.005'))**(m/math.mpf(2))*math.exp((-i*rs + v)/(ft*nr))/(ft*nr)

    fh = iss**XT*rs*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(math.exp((-i*rs + v)/(ft*n)) - XO)*math.exp((-i*rs + v)/(ft*n))/(XT*ft*n*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))) - iss*rs*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*math.exp((-i*rs + v)/(ft*n))/(ft*n)
    sh = isr*m*rs*(XO - (-i*rs + v)/vj)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(math.exp((-i*rs + v)/(ft*nr)) - XO)/(vj*((XO - (-i*rs + v)/vj)**XT + XF)) - isr*rs*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*math.exp((-i*rs + v)/(ft*nr))/(ft*nr)

    return np.matrix ([[fh+sh]])

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

    return np.matrix ([[iss*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(-math.exp((-i*rs + v)/(ft*n)) + XO)*(math.exp((-i*rs + v)/(ft*n)) - XO)/(XT*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))) + math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(math.exp((-i*rs + v)/(ft*n)) - XO),
                         iss**XT*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(-i*rs + v)*(math.exp((-i*rs + v)/(ft*n)) - XO)*math.exp((-i*rs + v)/(ft*n))/(XT*ft*n**XT*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))) - iss*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(-i*rs + v)*math.exp((-i*rs + v)/(ft*n))/(ft*n**XT),
                        iss*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))*(-ikf/(XT*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))**XT) + XO/(XT*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))))*(math.exp((-i*rs + v)/(ft*n)) - XO)/ikf,
                         ((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(math.exp((-i*rs + v)/(ft*nr)) - XO),
                         -isr*(-i*rs + v)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*math.exp((-i*rs + v)/(ft*nr))/(ft*nr**XT),
                        isr*m*(XO - (-i*rs + v)/vj)*(-i*rs + v)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(math.exp((-i*rs + v)/(ft*nr)) - XO)/(vj**XT*((XO - (-i*rs + v)/vj)**XT + XF)),
                         isr*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(math.exp((-i*rs + v)/(ft*nr)) - XO)*math.log((XO - (-i*rs + v)/vj)**XT + XF)/XT,
                         i*iss**XT*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*(math.exp((-i*rs + v)/(ft*n)) - XO)*math.exp((-i*rs + v)/(ft*n))/(XT*ft*n*(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO))) - i*iss*math.sqrt(ikf/(ikf + iss*(math.exp((-i*rs + v)/(ft*n)) - XO)))*math.exp((-i*rs + v)/(ft*n))/(ft*n)+ \
                         i*isr*m*(XO - (-i*rs + v)/vj)*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*(math.exp((-i*rs + v)/(ft*nr)) - XO)/(vj*((XO - (-i*rs + v)/vj)**XT + XF)) - i*isr*((XO - (-i*rs + v)/vj)**XT + XF)**(m/XT)*math.exp((-i*rs + v)/(ft*nr))/(ft*nr)
                        ]])

def Jac_Kirch_DiodeV2Mod2DirectBranchmathath (x,b,c,y):
    """
    [Реестровая]
    непосредственная проверка невозможна
    :param x:
    :param b:
    :param c:
    :param y:
    :return:
    """
    #return dfdy(y,x,b,c)**(-1) * dfdb(y,x,b,c)

    return np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


def solver_Kirch_DiodeV2Mod2DirectBranchmathath (x,b,c=None):
    """
    [Реестровая] +
    двухступенчатый - сначала np solver, затем math solver
    :param x:
    :param b:
    :param c:
    :param npVersion: если true, то возвращает numpy-compatible результат первой ступени, иначе идёт дальше
    :return:
    """
    global numnone
    global FT
    func = lambda y, x, b, c:  [float(func_Kirch_DiodeV2Mod2DirectBranchmathath (y,x,b,c)[0])]

    solvinit=[0.001]
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


    return solx






def test_Kirch_DiodeV2Mod2DirectBranchmathath():
    btrue=[341.4e-6, 2.664, 37.08e-3, 17.26e-27, 5.662, 4.282, 0.5751, 3.65e-3]

    btrueNormalDiode=[6.233e-8, 2.1, 1000, 17.26e-100, 2, 1, 0.5, 3.665e-2]

    #btrueNormalDiode=[math.mpf('7.69e-8'), math.mpf('1.45') , math.mpf('10'), math.mpf(0),math.mpf(2),math.mpf(1),math.mpf('0.5'), math.mpf('.0422')] #номинальные значения диода D1N4001 с сайта, вроде официальной модели производителя

    b=btrue
    rng=np.arange(0.001,1.5,0.01)
    Ve=np.array([ [1.9e-5] ]  )
    c=None

    #снимем ВАХ
    #resrng=[solver_Kirch_DiodeV2Mod2DirectBranchmathath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.

    resrng=list()
    nrng=list()

    resrng1,nrng1=list(), list()

    for x in rng:
        slv=solver_Kirch_DiodeV2Mod2DirectBranchmathath ([x],b)
        slv1=solver_Kirch_DiodeV2Mod2DirectBranchmathath ([x],btrueNormalDiode)

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

    #resrng=[o_pmath.makeMeasOneDot_lognorm(solver_Kirch_DiodeV2Mod2DirectBranchmathath, [x],b,c,Ve) for x in rng]


    #resrng1=[solver_Kirch_DiodeV2Mod2DirectBranchmathath ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.


    plt.plot(nrng1 , resrng1, label='normal')
    plt.plot(nrng , resrng, label='Schottky')
    plt.legend(loc='upper left')

    plt.axis([0,1,0,1])

    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()



def extraction_Kirch_DiodeV2Mod2DirectBranchmathath():
    """
    [Реестровая]

    :return:
    """
    global FT,XO,XT,XF
    global foldername

    funcf=solver_Kirch_DiodeV2Mod2DirectBranchmathath
    jacf = Jac_Kirch_DiodeV2Mod2DirectBranchmathath
    c=None
    Ve=np.array([ [1e-5] ]  )
    #       0           1       2        3         4       5    6        7
    btrue=[341.4e-6, 2.664, 37.08e-3, 17.26e-27, 5.662, 4.282, 0.5751, 3.65e-3]

    #коридор в 20 процентов в обе стороны
    diap = .1
    bstart = [item-item*diap for item in btrue]
    bend = [item+item*diap for item in btrue]

    binit=btrue

    xstart, xend = [0.001] , [1]

    N=20
    print("performing aprior plan:")

    import os
    filename ='N{0}_'.format(N)+os.path.basename(__file__).replace('.py','_plan')
    oplan = o_ap.makePlanCached (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6, verbose=True, foldername=foldername, cachname=filename, verbosePlan=True)


test_Kirch_DiodeV2Mod2DirectBranchmathath()
#extraction_Kirch_DiodeV2Mod2DirectBranchmathath()












