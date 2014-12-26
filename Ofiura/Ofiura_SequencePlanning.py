__author__ = 'vasilev_is'
__author__ = 'vasilev_is'
import copy
import math

import numpy as np

import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Estimation as o_e
import Ofiura.Ofiura_general as o_g
import Ofiura.Ofiura_ApriorPlanning as o_ap






#http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
def getbSeqPlanUltra (xstart:list, xend:list, N:int, btrue:list, binit:list, c, Ve, jacf, funcf, initplan=None, NSIG=10, smallestdetVb=1e-6, implicit=False, lognorm=False, dotlim=500, verbose=False, terminationOptDict={}):
    """
    Осуществляет последовательное планирование и оценку коэффициентов модели
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: количество точек в плане эксперимента
    :param binit: стартовое значение вектора коэффициентов
    :param btrue: истинное значение вектора коэффициентов
    :param c: словарь дополнительных переменных
    :param Ve:  ковариационная матрица y
    :param jacf: функция, возвращающая якобиан модели, входящие x,b,c,y
    :param funcf: функция, возвращающая значение модели, входящие x,b,c
    :param initplan: - начальный план
    :param NSIG:
    :param implicit - неявность функции. Если функция неявна, то в gknu sign=0
    :param lognorm - требуется ли использование логнорм для получения измеренных данных
    :param dotlim - предел добавленных точек
    :param verbose - выводить ли информацию по итерациям
    :param terminationOptDict - словарь опций завершения итеративного процесса:
        VdShelfPow - по умолчанию -1 math.fabs(detVb-prevdetVb)/prevdetVb<math.pow(10,-VdShelfPow) NSIG для полки по detVb

    :return: k, число итераций, лог
    Для переделывания на реальные измерения btrue=None и функции моделирования переводятся на измерительные // btrue вообще должен бы быть исключен из всех функций ofiura
    """


    startplan = initplan if initplan else o_p.makeRandomUniformExpPlan(xstart, xend, N)
    origplan = copy.copy(startplan)
    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c, Ve) if lognorm else o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    log=""
    b=binit
    prevdetVb=None
    for numiter in range(dotlim): #ограничитель цикла - если выход произошёл по ограничению, значит, возможна ошибка
        estim=o_e.grandCountGN_UltraX1(funcf, jacf, measdata, b, c, NSIG, implicit=implicit, verbose=False) #получили оценку b binit=b
        #measdata почему-то разная, причины неизвестны.
        b=estim[0]
        Sk=estim[1]
        Vb=o_p.countVbForMeasdata(b,  c, Ve, jacf, measdata)
        #посчитали определитель
        detVb=np.linalg.det(Vb)

        if verbose:
            print ("Sequence Plan Iteration: {0}\nb={1}\ndetVb={2}\nprevdetVb={3} \nSk={4}".format(numiter, b, detVb, prevdetVb, Sk))


        VdShelfPow=terminationOptDict['VdShelfPow'] if 'VdShelfPow' in terminationOptDict else -1

        condition=prevdetVb!=None and math.fabs(detVb-prevdetVb)/prevdetVb<math.pow(10,VdShelfPow) #если вышли на плато



        prevdetVb=detVb

        if condition:
            return b, numiter, Sk, startplan, origplan, log, estim, measdata, detVb #то всё вернуть

        #иначе поиск следующей точки плана

        xdot=copy.copy(xstart) #получили начальную точку начальное значение - ровно пополам диапазона
        for i in range(len(xstart)): #присвоили ей значение
            xdot[i]=xstart[i]+(xend[i]-xstart[i])/2

        #объектная функция
        measure=lambda x,b,c:o_p.makeMeasOneDot_lognorm(funcf, x, b, c, Ve) if lognorm else o_p.makeMeasOneDot(funcf, xdot, b, c, Ve)

        function = lambda x: np.linalg.det(o_p.countVbForPlan(o_g.appendToList(startplan, x),b,c,Ve,jacf, measure))
        #function = lambda x: np.linalg.det(o_p.countVbForPlan(o_g.appendToList(startplan, x),b,c,Ve,jacf, funcf))
        #funcf заменено на measure, которая добавляет дисперсию.

        #function и measure есть разные функции - первая даёт идеальный результат, вторая - с дисперсией
        #создать функцию, которая будет возвращать полученные от func данные, налагая дисперсию.
        #каждый раз будет пытаться добавить в план точку и вернуть определитель с добавленной точкой
        #где тут добавлять дисперсию, а где - не добавлять, есть тащем-та вопрос открытый

        xdot=o_g.doublesearch (xstart, xend, xdot, function) #оптимизировали значение точки

        startplan.append(xdot)
        #measdata.append({'x':xdot, 'y':funcf(xdot,b,c)})

        ymeas=o_p.makeMeasOneDot_lognorm(funcf, xdot, btrue, c, Ve) if lognorm else o_p.makeMeasOneDot(funcf, xdot, btrue, c, Ve)

        #примитивная защита от кривых результатов измерения
        if ymeas!=None:
            measdata.append({'x':xdot, 'y':ymeas})

        #measdata = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c, Ve) if lognorm else o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)


    #окончание этого цикла "естественным путём" говорит о том, что превышено максимальное число итераций
    return b, dotlim, Sk, startplan, origplan, log+"ERROR: maximum number of iterations archieved", estim, measdata, detVb


def getHybridPlan (xstart:list, xend:list, N:int, btrue:list, binit:list, bstart:list, bend:list, c, Ve, jacf, funcf, NSIG=10, implicit=False, lognorm=False, dotlim=500, verbose=False):
    oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, func=None, Ntries=5)[1] #сперва строим априорный план
    return getbSeqPlanUltra (xstart, xend, N, btrue, binit, c, Ve, jacf, funcf, initplan=oplan, dotlim=500) #теперь последовательный с лимитом на 500 точек лишних



