__author__ = 'vasilev_is'

import copy
import random
import math

import numpy as np
from scipy.optimize import minimize

"""
replaceInList (list, i, x):
    :param list: список
    :param i: индекс элемента
    :param x: переменная, на которую заменяем элемент
    :return: Возвращает список list, в котором i-й элемент заменяется на x
 uniformVector(xstart, xend):
    :param xstart: вектор начала диапазона
    :param xend: вектор конца диапазона
    :return: вектор, являющийся частью равномерного распределения между xstart и xend
makeUniformExpPlan(xstart:list, xend:list, N:int):
    Делает равномерный неслучайный план эксперимента
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param N: Количество точек в плане
    :return: равномерный план эксперимента

    makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    Делает случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param N: Количество точек в плане


     makeMeasAccToPlan(func, expplan:list, b:list, c:dict, ydisps:list=None, n=1, outfilename="", listOfOutvars=None):
    Моделирует измерения в соответствии с планом эксперимента.
    :param func: векторная функция, на вход принимает x, на выходе выдаёт значения y
    :param expplan: план эксперимента (список значений вектора x) (список списков)
    :param b: вектор b - коэффициенты
    :param c: словарь c - дополнительные переменные
    :param ydisps: вектор дисперсий y (диагональ ковариационной матрицы)
    :param n: объём выборки y  - в данной версии не используется
    :param outfilename: имя выходного файла, куда писать план - в данной версии не используется
    :param listOfOutvars: список выносимых переменных - в данной версии не используется
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...

    def countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    Считает значение определителя Vb, ковариационной матрицы коэффициентов, для данного плана эксперимента
    :param expplan: план эксперимента (список списков)
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента

    def countMeanVbForAprior_S4000(expplan:list, bstart:list, bend:list, c, Ve, jac, func=None):
    Считает среднее значение определителя и дисперсию определителя Vb, формируя выборку из значений случайных b, распределённых равномерно (uniform) в диапазоне. Выборка в данной версии объёмом в 30 значений.
    :param expplan: план эксперимента
    :param bstart: начало диапазона b (вектор)
    :param bend: конец диапазона b (вектор)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов
    :param jac: функция якобиана (на входе x,b,c=None, y=None)
    :return: среднее значение определителя [0] и его дисперсию [1]

    doublesearch (xstart, xend, xinit, function):
    Реализует метод двойного поиска с ограничениями. Внимание: функция требует правки в части штрафов в зависимости от значений, которые может принимает объектная функция.
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param xinit:начальное значение x (вектор)
    :param function: объектная функция, принимает на вход x
    :return: оптимизированное значение x

    def grandApriornPlanning (xstart:list, xend:list, N:int, bstart:list, bend:list, c, Ve, jac, func=None, Ntries=30):
    Реализует априорное планирование эксперимента
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x  (вектор)
    :param N: размер плана (количество контрольных точек)
    :param bstart: начало диапазона b (вектор)
    :param bend: конец диапазона b (вектор)
    :param c: словарь дополнительных переменных
    :param Ve: Ковариационная матрица y, реально её диагональ (вектор)
    :param jac: Якобиан функции, принимает на вход x,b,c,y
    :return: кортеж: 0: оптимизированное значение определителя Vb, 1: оптимальный план эксперимента

"""

def replaceInList (list, i, x):
    """
    :param list: список
    :param i: индекс элемента
    :param x: переменная, на которую заменяем элемент
    :return: Возвращает список list, в котором i-й элемент заменяется на x
    """
    l=copy.deepcopy(list)
    l[i]=x
    return l   #+

def uniformVector(xstart, xend):
    """
    :param xstart: вектор начала диапазона
    :param xend: вектор конца диапазона
    :return: вектор, являющийся частью равномерного распределения между xstart и xend
    """
    res= [0 for i in range(len(xstart))]
    for i in range(0, len(xstart)):
        res[i]= random.uniform (xstart[i], xend[i])
    return res
    #Можно сделать map для списка списков сразу (или просто для нескольких списков) на основе этой функции  #

def makeUniformExpPlan(xstart:list, xend:list, N:int):
    """
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: равномерный план эксперимента
    """

    res=list()

    xstartnp=np.array(xstart)
    xendnp=np.array(xend)

    xstep = (xendnp-xstartnp)/N

    for i in range(N):
        res.append(list(xstart+i*xstep))

    return res

def makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    """
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор
    """
    res=list()
    for i in range(0, N):
        res.append(uniformVector(xstart, xend))
    return res

def makeMeasAccToPlan(func, expplan:list, b:list, c:dict, ydisps:list=None, n=1, outfilename="", listOfOutvars=None):
    """

    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param ydisps: вектор дисперсий y (диагональ ковариационной матрицы)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """


    res = list()

    for i in range(len(res)):
        y=func(expplan[i],b,c)
        #Внесём возмущения:
        if not ydisps==None:
            for k in range(len(y)):
                y[k]=random.normalvariate(y[k], math.sqrt(ydisps[k]))


        curdict = {'x':expplan[i], 'y':y}
        #res[i]["y"]=y
        res.append(curdict)
    return res

def countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    """
    :param expplan: план эксперимента
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента
    """
    G=np.zeros((len(b),len(b))) #матрица G

    for point in expplan:
        jj=jac(point, b, c, func(point,b,c) if func else None)
        #G+=jj*np.linalg.inv(Ve)*jj.T
        G+=jj*jj.T
        #print (jj*jj.T)

    return np.linalg.inv(G)

def countMeanVbForAprior_S4000(expplan:list, bstart:list, bend:list, c, Ve, jac, func=None):
    """

    :param expplan: план эксперимента
    :param bstart: начало диапазона b (вектор)
    :param bend: конец диапазона b (вектор)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов
    :param jac: функция якобиана (на входе x,b,c=None, y=None)
    :return: среднее значение определителя [0] и его дисперсию [1]
    """
    DS=0 #среднее определителя
    SD=0 #дисперсия определителя

    for sss in range(1, 30): #30 - количество  попыток в выборке
        b=uniformVector (bstart, bend)
        Vb=countVbForPlan(expplan, b, c, Ve, jac, func)
        D=np.linalg.det(Vb)


        if D:
            DS=(D+(sss-1)*DS)/sss  #среднее определителя
            SD=((DS-D)*(DS-D)+(sss-1)*SD)/sss #дисперсия определителя


    return DS, SD

def doublesearch (xstart, xend, xinit, function):
    """
    Метод двойного поиска с ограничениями
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param xinit:начальное значение x (вектор)
    :param function: объектная функция, принимает на вход x
    :return: оптимизированное значение x
    """
    x=copy.deepcopy(xinit)
    #xcurr=copy.deepcopy(x)
    for i in range (len(x)):
        step=(xend[i]-xstart[i])/3 #задание начального шага
        brval=10000000 #критичное число итераций, после которого надо выйти из цикла и сообщить об ошибке
        while(1):
            xcurr=copy.deepcopy(x)
            d3=function(xcurr)
            xcurr[i]=x[i]+step  #пробуем шаг вперёд
            d4=function(xcurr) if xcurr[i]<xend[i] else 1 #штраф при выходе из диапазона
            xcurr[i]=x[i]-step #пробуем шаг назад
            d5=function(xcurr) if xcurr[i]>xstart[i] else 1 #штраф при выходе из диапазона

            # print ('head')
            # print (d3,d4,d5, step, x, i)
            # print ('tail')


            if (d3-d4)*(d3-d5)>0: #
                step=step/2       #уменьшение шага вполовину
                if step<(xend[i]-xstart[i])/1000: #если шаг меньше 1/1000, значит этот член вектора прооптимизировали, переходим к следующему
                    break
            else:
                if (d4>d5):
                    x[i]=x[i]-step #шаг назад
                else:
                    x[i]=x[i]+step #шаг вперёд

            brval-=1
            if not brval:
                print ("doublesearch: exit due to max. number of iteration archieved, possible ERROR")
                break

    return x

def grandApriornPlanning (xstart:list, xend:list, N:int, bstart:list, bend:list, c, Ve, jac, func=None, Ntries=30):
    """
    Реализует априорное планирование эксперимента
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x  (вектор)
    :param N: размер плана (количество контрольных точек)
    :param bstart: начало диапазона b (вектор)
    :param bend: конец диапазона b (вектор)
    :param c: словарь дополнительных переменных
    :param Ve: Ковариационная матрица y, реально её диагональ (вектор)
    :param jac: Якобиан функции, принимает на вход x,b,c,y
    :return: кортеж: 0: оптимизированное значение определителя Vb, 1: оптимальный план эксперимента
    """

    dopt=100000000
    planopt=None

    for i in range(0,Ntries):
        plan = makeRandomUniformExpPlan(xstart, xend, N)


        unopt=countMeanVbForAprior_S4000(plan, bstart, bend, c, Ve, jac)[0]




        #оптимизация
        for j in range(N):
            xdot=copy.deepcopy(plan[j])
            function = lambda x: countMeanVbForAprior_S4000(replaceInList(plan,j,x), bstart, bend, c, Ve, jac, func)[0]   #TODO test it!!

            # boundsarr=list()
            # for k in range(len(xstart)):
            #     boundsarr.append((xstart[k],xend[k]))
            #
            # sol = minimize (function, xdot, bounds=boundsarr)
            # plan[j]=sol.x
            plan[j]=doublesearch(xstart, xend, xdot, function)



        dcurr=countMeanVbForAprior_S4000(plan, bstart, bend, c, Ve, jac)[0]

        print ("unoptimized-optimized:",unopt, dcurr)

        if (dcurr<dopt):
            dopt=dcurr
            planopt=plan


    return (dopt, planopt)






def test():
    """
    Тестирует априорное планирование
    :return:
    """
    xstart=[1, 100]
    xend=[20,200]
    N=10
    c={"a":1000}
    funcf=lambda x,b,c: np.array ( [ b[1]+b[2]*x[1]+b[3]*x[2]+b[4]*x[1]*x[2]+b[5]*x[1]*x[1]+b[6]*x[2]*x[2],   b[7]+b[8]*x[1]+b[9]*x[2]+b[10]*x[1]*x[2]+b[11]*x[1]*x[1]+b[12]*x[2]*x[2] ] )
    jacf = lambda x,b,c,y: np.matrix([ [1, x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1], 0, 0, 0, 0, 0, 0],
                                       [0,0,0,0,0,0,1,x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1]] ]).T
    Ve=np.array([ [0.1, 0],
                  [0, 0.1]]  )
    bstart=[0.8,0.4,1.4,0.2,0.9,0.3,1.4,0.2,2.1,3.1,4.1,5.1]
    blen=  [0.3,0.2,0.2,0.2,0.2,0.3,0.2,0.2]
    bend=  [1.1,0.6,1.6,0.4,1.1,0.6,1.6,0.4,2.5,3.3,4.6,5.6]

    #print (doublesearch ([1, 0.5], [10,10], [9,9], lambda x: x[0]*x[0]+2*x[1]*x[1]+10)) #тестирование поиска

    rs=grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, func=None, Ntries=30)
    print (rs[0])

    print ('Experimental plan')
    for r in rs[1]:
        print(r[0], '\t', r[1])


    # plan=makeUniformExpPlan(xstart, xend, N)
    # func = lambda x,b,c: [x[0]*b[0]+c["a"], x[1]*b[1]+c["a"], x[2]*b[2]+c["a"]]
    # meas = makeMeasAccToPlan(func, plan,  b, c, [0.0001]*3)
    # for x in meas:
    #     print (x)

