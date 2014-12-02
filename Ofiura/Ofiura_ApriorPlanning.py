from Ofiura import Ofiura_planning as o_p

__author__ = 'vasilev_is'
import copy

import numpy as np

import Ofiura.Ofiura_general as o_g
import traceback

"""
def countMeanVbForAprior_S4000(expplan:list, bstart:list, bend:list, c, Ve, jac, func=None):
    Считает среднее значение определителя и дисперсию определителя Vb, формируя выборку из значений случайных b, распределённых равномерно (uniform) в диапазоне. Выборка в данной версии объёмом в 30 значений.
    :param expplan: план эксперимента
    :param bstart: начало диапазона b (вектор)
    :param bend: конец диапазона b (вектор)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов
    :param jac: функция якобиана (на входе x,b,c=None, y=None)
    :return: среднее значение определителя [0] и его дисперсию [1]


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
        b=o_g.uniformVector (bstart, bend)
        Vb=o_p.countVbForPlan(expplan, b, c, Ve, jac, func)
        D=np.linalg.det(Vb)


        if D:
            DS=(D+(sss-1)*DS)/sss  #среднее определителя
            SD=((DS-D)*(DS-D)+(sss-1)*SD)/sss #дисперсия определителя


    return DS, SD

def grandApriornPlanning (xstart:list, xend:list, N:int, bstart:list, bend:list, c, Ve, jac, func=None, Ntries=30, verbosePlan=False, initplan=None, verbose=False):
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
    :param verbosePlan: если true, пишет все планы в файлы, иначе только оптимальный
    :param initplan: начальный план. По умолчанию случаен, но может быть задан (в этом случае делается одна попытка)
    :param verbose: выдавать ли в консоль информацию optimized-original
    :return: кортеж: 0: оптимизированное значение определителя Vb, 1: оптимальный план эксперимента
    """

#Апдейт: теперь по умолчанию в список планов заряжается равномерный

    dopt=100000000
    planopt=None

    Ntries1=Ntries+1
    if initplan!=None:
        Ntries1=1


    if verbose:
        print('\n\nДанные априорного планирования:')
        print('Неоптимизированное-оптимизированное значение среднего det(Vb)')

    for i in range(0,Ntries1):
        try:
            if initplan==None:
                m=len(xstart) #длина вектора входных параметров
                plan = o_p.makeUniformExpPlan(xstart, xend, N**(1/float(m))) if i==0 else o_p.makeRandomUniformExpPlan(xstart, xend, N)

                if verbose:
                    print('plan length:', len(plan))

            else:
                plan = initplan
            unopt=countMeanVbForAprior_S4000(plan, bstart, bend, c, Ve, jac, func)[0]
            #оптимизация
            for j in range(N):
                xdot=copy.deepcopy(plan[j])
                function = lambda x: countMeanVbForAprior_S4000(o_g.replaceInList(plan,j,x), bstart, bend, c, Ve, jac, func)[0]

                boundsarr=list()
                for k in range(len(xstart)):
                    boundsarr.append((xstart[k],xend[k]))

                #sol = optimize.minimize (function, xdot, bounds=boundsarr)
                #plan[j]=sol.x
                #В этом варианте не работает, хоть результат и проходит быстрее

                plan[j]=o_g.doublesearch(xstart, xend, xdot, function)



            dcurr=countMeanVbForAprior_S4000(plan, bstart, bend, c, Ve, jac, func)[0]

            if verbose:
                print ("{0} unoptimized-optimized:".format('uniform' if i==0 else ''),unopt, dcurr)
            if dcurr<dopt or planopt==None:
                dopt=dcurr
                planopt=plan


            if (verbosePlan):
                o_p.writePlanToFile(plan, "{0}plan.txt".format(i))

        except BaseException as e:
            print ('This try failed, due to exception e=',e)
            tb = traceback.format_exc()
            print(tb)

    return dopt, planopt

