__author__ = 'vasilev_is'

import random
import math

import numpy as np

import Ofiura_general as o_g


"""

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

countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    Считает значение определителя Vb, ковариационной матрицы коэффициентов, для данного плана эксперимента
    :param expplan: план эксперимента (список списков)
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента

writePlanToFile (plan, filename='plan.txt'):
    Пишет план в файл. Значения пишутся в формате python (вектор [])
    :param plan: план (список)
    :param filename: имя файла
    :return: ничего

readPlanFromFile (filename='plan.txt'):
    Читает план из файла - запись в виде столбца векторов []
    :param filename: имя файла
    :return: план (список векторов)

def makeUniformExpPlan(xstart:list, xend:list, N:int):
    Создаёт равномерный план эксперимента
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане (внимание! это количество измерений каждого вектора, то есть, реальное кол-во будет N^len(xstart))
    :return: равномерный план эксперимента

makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор

makeMeasAccToPlan(func, expplan:list, b:list, c:dict, Ve=[], n=1, outfilename="", listOfOutvars=None):
    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...


countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    :param expplan: план эксперимента
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента

"""


def for_filter (x,lim):
    for val in x['y']:
        if val>lim:
            return False
    return True

def filterList(inmeasdata, lim=1e55):
    filt = lambda x: for_filter(x,lim)
    return list(filter(filt, inmeasdata))






def writePlanToFile (plan, filename='plan.txt'):
    """
     Пишет план в файл. Значения пишутся в формате python (вектор [])
    :param plan: план (список)
    :param filename: имя файла
    :return: ничего
    """
    with open (filename, 'wt') as outfile:
        for point in plan:
            outfile.write(point.__str__())
            outfile.write('\n')

def readPlanFromFile (filename='plan.txt'):
    """
    Читает план из файла - запись в виде столбца векторов []
    :param filename: имя файла
    :return: план (список векторов)
    """
    with open (filename, 'rt') as infile:
        lines = infile.readlines()
        return list(map(eval, lines))

def makeUniformExpPlan(xstart:list, xend:list, N:int):
    """
    Создаёт равномерный план эксперимента
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане (внимание! это количество измерений каждого вектора, то есть, реальное кол-во будет N^len(xstart))
    :return: равномерный план эксперимента
    """
    res=list()

    xstartnp=np.array(xstart)
    xendnp=np.array(xend)
    xstep = list((xendnp-xstartnp)/N)

    evalstr="import numpy\n\n"
    lststr=""

    for i in range (len(xstart)):
        evalstr+="\t"*i+"for x{0} in numpy.arange(xstart[{0}], xend[{0}], xstep[{0}]):\n".format(i)
        lststr+="x{0}".format(i)+("," if i+1<len(xstart) else "")
    evalstr+="\t"*(i+1)+"res.append(("+lststr+"))"

    #print (evalstr)
    exec(evalstr,  locals()) #исполняет полученную сроку, собсна, формирует список входных переменных

    # for i in range(N):
    #     res.append(list(xstart+i*xstep))
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
        res.append(o_g.uniformVector(xstart, xend))
    return res

def makeMeasAccToPlan(func, expplan:list, b:list, c:dict, Ve=[], n=1, outfilename="", listOfOutvars=None):
    """

    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """
    res = list()

    for i in range(len(expplan)):
        y=func(expplan[i],b,c)
        #Внесём возмущения:
        if not Ve==None:
            ydisps=np.diag(Ve)
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
        jj=np.array(jac(point, b, c, func(point,b,c) if func else None))
        #G+=jj*np.linalg.inv(Ve)*jj.T

        G+=np.dot(jj.T, jj)


    return np.linalg.inv(G)


