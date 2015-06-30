__author__ = 'reiner'

import ModelMaker.MMModel as mmd
import ModelMaker.MMMeasdata as mmsd

import Ofiura.Ofiura_planning as op
import Ofiura.Ofiura_Estimation as o_e

"""
Главный управляющий файл
+++  Тестирование  +++
Общий алгоритм:
1. Родить регрессионную модель
2. Установить btrue и Ve
3. Сделать measdata
4. Оценить коэффициенты регресии
"""

def test ():
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






test()




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

