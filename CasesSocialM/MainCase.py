__author__ = 'reiner'
import scipy.optimize as scp

import ModelMaker.MMModel as mmd
import ModelMaker.MMMeasdata as mmsd

import Ofiura.Ofiura_planning as op
import Ofiura.Ofiura_Estimation as o_e
from statistics import mean


import numpy as np
import copy
import math

"""
Главный управляющий файл
+++  Тестирование  +++
Общий алгоритм:
1. Родить регрессионную модель
2. Установить btrue и Ve
3. Сделать measdata
4. Оценить коэффициенты регресии
"""

def mainModellingPrc ():
    measdata = mmsd.MMMeasdata('Table.csv')
    pool=['In1.1', 'In1.2', 'In1.3', 'In1.4', 'In2.1', 'In2.2', 'In2.3', 'In5.1', 'In5.2', 'In5.3', 'In5.4']
    modelcls=mmd.MMModel()
    modelcls.makeLinearRegression(len(pool))
    binit=[1]*(len(pool)+1)

    #measdata = op.makeMeasAccToPlan_lognorm(modelcls.solver,plan,btrue,None,Ve)

    measdata.inlist=pool
    measdata.outlist=['O1.1']
    measdata.showReadableNames(1,1)

    rs=o_e.grandCountGN_UltraX1(modelcls.solver,modelcls.jacf,measdata,binit,None,NSIG=10,implicit=1)
    print (rs)


    averages= [134676.99,	20.63,	2.57,	11.44,	25.78,	128,	1134.7,	-9.76,	9.29,	64,	0.315,	11.26,	5.638,	41.6,	34.787	,109.32	,50	,49.4]



def ifgrp (instr:str, grp:int, inoit='In'):
    instr1 = instr.replace('In','').replace('O','')
    return (int(instr1.split('.')[0])==grp) and (inoit in instr)

def util():
    #описание пулов исходных данных по подгруппам
    #рассматриваем только группу "Благосостояние"
    ps=[[1,2,5],    #1
        [1,3,5],    #2
        [1],        #3
        [2,6],      #4
        [6],        #5
        [1,3,6],    #6
        [1,3,6]]    #7

    measdata = mmsd.MMMeasdata('Table.csv')

    print (measdata.ids)

    for pp in ps:
        inpool=[]
        for i in pp:
            inpool+=[item for item in measdata.ids if ifgrp(item,i)]

        print (inpool)

#исходные пулы для переменных
# ['In1.1', 'In1.2', 'In1.3', 'In1.4', 'In2.1', 'In2.2', 'In2.3', 'In5.1', 'In5.2', 'In5.3', 'In5.4']
# ['In1.1', 'In1.2', 'In1.3', 'In1.4', 'In3.1', 'In3.2', 'In3.3', 'In3.4', 'In3.5', 'In5.1', 'In5.2', 'In5.3', 'In5.4']
# ['In1.1', 'In1.2', 'In1.3', 'In1.4']
# ['In2.1', 'In2.2', 'In2.3', 'In6.1', 'In6.2', 'In6.3', 'In6.4', 'In6.5', 'In6.6', 'In6.7']
# ['In6.1', 'In6.2', 'In6.3', 'In6.4', 'In6.5', 'In6.6', 'In6.7']
# ['In1.1', 'In1.2', 'In1.3', 'In1.4', 'In3.1', 'In3.2', 'In3.3', 'In3.4', 'In3.5', 'In6.1', 'In6.2', 'In6.3', 'In6.4', 'In6.5', 'In6.6', 'In6.7']
# ['In1.1', 'In1.2', 'In1.3', 'In1.4', 'In3.1', 'In3.2', 'In3.3', 'In3.4', 'In3.5', 'In6.1', 'In6.2', 'In6.3', 'In6.4', 'In6.5', 'In6.6', 'In6.7']






#
    # fr=measdata[0]
    #
    # with open('test1.csv','wt') as csvf:
    #     lst = ['x{0}'.format(i) for i in range(len(fr['x']))]
    #     lst+=[ 'y{0}'.format (j) for j in range(len(fr['y']))]
    #     csvf.write(','.join(lst))
    #     csvf.write('\n')
    #     #читаемые идентификаторы вбиваются ручками
    #
    #     for point in measdata:
    #         jlst = point['x'] + point['y']
    #         line=','.join([str(h) for h in jlst])
    #         csvf.write(line)
    #         csvf.write('\n')


def estimateLinRegrPool (measdata):
    """
    функция, принимает на вход measdata, строит модель по длине входных переменных,
    определяет коэффициенты линейной регрессии
    :param measdata:
    :param pool:
    :return:
    """
    pool=measdata.inlist

    modelcls=mmd.MMModel()
    modelcls.makeLinearRegression(nvar=len(pool))
    binit=[1]*(len(pool)+1)
    xarr=np.matrix(measdata.getXarray()).T
    yarr=measdata.getY()
    bestinit=scp.curve_fit(modelcls.solverSc,xarr,yarr,p0=binit)[0]
    Skbestinit = sum (   [ (modelcls.solver(point['x'],bestinit)[0]-point['y'][0])**2 for point in measdata ])
    return bestinit, Skbestinit







def simpleTest(measdata, pool):
    """
    отбирает наиболее важные факторы из пула, строит линейную регрессию, выводит все данные.
    Первый компонент вектора регрессии - свободный член.
    :param measdata:
    :param pool:
    :return:
    """

    measdata.inlist=pool
    bestinit, Skbestinit= estimateLinRegrPool(measdata)
    print ('Initial guess and Sk')
    print (bestinit, Skbestinit)


    minipool=copy.copy(pool)

    varskdict={}
    for var in pool:
        minipool.remove(var)
        print (var)

        measdata.inlist=minipool
        best, Skbest= estimateLinRegrPool(measdata)

        print (best, Skbest)

        varskdict[var]=Skbest
        minipool.append(var)

    print ('SK for pool without this variable')
    print (varskdict)
    print ('Diffence with initial Sk, means importance of the variable')
    importance  = {}
    for k, v in varskdict.items():
        importance[k]=math.fabs(v-Skbestinit)
    print (importance)


    for var,sk in varskdict.items(): #перебираем все переменные
        if math.fabs(Skbestinit-sk)<3: #типа критическое число статистики хиквадрат уровень 0,05 степени свободы 9, по факту может надо степень свободы 1, тогда 0,0039
            minipool.remove(var) #если меньше, то вытряхнуть данную переменную

    print (minipool)
    #распечатываем пул переменных, который получился после фильтрования

    print ('estimation via reduced model')

    measdata.inlist=minipool
    print (estimateLinRegrPool(measdata))




measdata = mmsd.MMMeasdata('test1.csv')
pool=['in0','in1','in2','in3','in4','in5']
measdata.outlist=['out0']

simpleTest(measdata,pool)





def factorSelectionTest():

    measdata = mmsd.MMMeasdata('testData.csv')
    pool=['in1','in2','in3','in4','in5']
    measdata.outlist=['o1']

    modelcls=mmd.MMModel()

    #1 посчитать среднее по выходной переменной
    ylist = measdata.getCut('o1')
    print (ylist)
    avy=mean(ylist)
    #это и является самой примитивной, вырожденной регрессией

    Skinit=sum([(y-avy)**2 for y in ylist])
    print ('Initial Sk: {0}'.format(Skinit))

    for nvars in range(1,len(pool)):
        #для числа переменных от 1 до длины пула
        modelcls.makeLinearRegression(nvars) #построим регрессию с 1 переменной
        binit=[1]*(nvars+1)



        SkCurr=Skinit
        AddingVarcurr=''
        rscurr=None

        #measdata.inlist.append(pool[0]) #сразу всунем первую переменную, чтоб потом вытягивать её, если шо

        for var in pool: #для каждой переменной в пуле
            if measdata.inlist:
                measdata.inlist.pop()
            measdata.inlist.append(var) #заменить последнюю переменную такой

            #вывод отладочных данных
            print ('Model, variables number={0}'.format(nvars))
            measdata.showReadableNames(1,1)

            rs1=o_e.grandCountGN_UltraX1(modelcls.solver,modelcls.jacf,measdata,binit,None,NSIG=5,implicit=1, maxiter=10)
            rs2=o_e.grandCountGN_UltraX1(modelcls.solver,modelcls.jacf,measdata,binit,None,NSIG=5,implicit=0, maxiter=10)
            #здесь присутствует некоторая проблема - implicit может и на explicit поменяться, фиг его
            # объектная (Sk и есть наше MCLL - Maximizing Component of Logarithm of Likelihood function)
            if rs1[4]>rs2[4]: rs=rs2
            else: rs=rs1

            if rs[4]<SkCurr: #если полученная на шаге Sk меньше, чем текущая. Здесь разница должна сравниваться по статистике! В таком исполнении последовательность переменных будет по степени важности
                print ('difference - LR statistics value {0}'.format(SkCurr-rs[4]))
                SkCurr=rs[4]
                AddingVarcurr=var
                rscurr=rs

        measdata.inlist.pop()

        if AddingVarcurr: #если есть, что добавлять
            print ('added a variable in regression')
            measdata.inlist.append(AddingVarcurr)
            pool.remove(AddingVarcurr)
        else:
            print ('finished, no significant variables left')
            print (rscurr)



    print (rs)


