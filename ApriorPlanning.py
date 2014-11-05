__author__ = 'vasilev_is'



import numpy as np

import Ofiura_ApriorPlanning as o_ap

import Ofiura_planning as o_p
import Ofiura_Estimation as o_e
import Ofiura_Qualitat as o_q


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
                                       [0,0,0,0,0,0,1,x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1]] ])

    #убрал .T

    Ve=np.array([ [0.1, 0],
                  [0, 0.1]]  )
    bstart=[0.8,0.4,1.4,0.2,0.9,0.3,1.4,0.2,2.1,3.1,4.1,5.1]
    blen=  [0.3,0.2,0.2,0.2,0.2,0.3,0.2,0.2]
    bend=  [1.1,0.6,1.6,0.4,1.1,0.6,1.6,0.4,2.5,3.3,4.6,5.6]

    #print (doublesearch ([1, 0.5], [10,10], [9,9], lambda x: x[0]*x[0]+2*x[1]*x[1]+10)) #тестирование поиска

    rs=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, func=None, Ntries=1)
    print (rs[0])

    print ('Experimental plan')
    for r in rs[1]:
        print(r[0], '\t', r[1])

    print(rs[1])

    measdata = o_p.makeMeasAccToPlan(funcf, rs, bend, c, Ve)
    #надо добавить скажем априорный план, с фильтрованием точек

    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, bstart, c,NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))

    print (gknu[0])




    # plan=makeUniformExpPlan(xstart, xend, N)
    # func = lambda x,b,c: [x[0]*b[0]+c["a"], x[1]*b[1]+c["a"], x[2]*b[2]+c["a"]]
    # meas = makeMeasAccToPlan(func, plan,  b, c, [0.0001]*3)
    # for x in meas:
    #     print (x)
test()




