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
    modelcls=mmd.MMModel()
    modelcls.makeLinearRegression(3)

    plan = op.makeUniformExpPlan([1,1,1],[4,4,4],3)

    btrue=[3,50,100,6]
    binit=[1,1,1,1]

    Ve=[[.1e-5]]
    
    #measdata = op.makeMeasAccToPlan_lognorm(modelcls.solver,plan,btrue,None,Ve)

    measdata = mmsd.MMMeasdata('test1.csv')
    measdata.inlist=['x0', 'x1', 'x2']
    measdata.outlist=['y0']


    measdata.showReadableNames(1)

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





    rs=o_e.grandCountGN_UltraX1(modelcls.solver,modelcls.jacf,measdata,binit,None,implicit=False)

    print (rs)






test()
