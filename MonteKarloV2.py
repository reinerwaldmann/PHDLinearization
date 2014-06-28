__author__ = 'vasilev_is'

import numpy as np

import sequence_generation as sg


#Monte-Karlo method applied to the resistor problem, using new sequence generator




def MonteKarloGenericV2 (funcstrdict, xvectorlistsdict, spreadvarslist, Vx, nvolx):
    """
    Определяет средние и дисперсию методом Монте-Карло
    """

    resdict=sg.generateAsDict (funcstrdict, xvectorlistsdict, spreadvarslist, Vx, nvolx, nvoly=1)
    meandict=dict() #средние
    vardict=dict() #дисперсия

    for key in resdict.keys():
        meandict[key]=np.mean(resdict[key])
        vardict[key]=np.var(resdict[key])

    print ("Mean values of outer variables:")
    for key in funcstrdict.keys():
        print (key, "\t", meandict[key])
    print ("VAR of outer variables:")
    for key in funcstrdict.keys():
        print (key, "\t", vardict[key])




#test area


funcstrdict= {"y1":"u1* (r2+r3)/(r1+r2+r3)", "y2":"u1* r3/(r1+r2+r3)"}
xvectorlistsdict = {"u1":[100],  "r1":[20], "r2":[30], "r3":[400]}
spreadvarslist  = ["r1", "r2", "r3"]
V=np.array       ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])
MonteKarloGenericV2 (funcstrdict, xvectorlistsdict, spreadvarslist, V, nvolx=10000)






