
__author__ = 'reiner'
import matplotlib.pyplot as plt
import numpy as np

import sequence_generation as sg

import ex3delitel as ex3




def mntkrleq():
    #в предположении,что имеет смысл держать центр эллипса совмещённым с центром прямоугольника
    global funcstrdict

    global spreadvarslist, V
    ksu2=0.955
    ksu3=0.888

    diap=0.001
    leftksu2=ksu2-diap
    rightksu2=ksu2+diap
    leftksu3=ksu3-diap
    rightksu3=ksu3+diap

    bestpercent=0
    vectr=[]

    for r2 in range(1,500, 1):
        r3=ksu3*r2/(ksu2-ksu3)
        r1=r2/(ksu2-ksu3)-r2-r3
        xvectorlistsdict = {"u1":[100],  "r1":[r1], "r2":[r2], "r3":[r3]}
        xvectorlistsdictc = {"u1":100,  "r1":r1, "r2":r2, "r3":r3}
        #для проверки
        centerksu2=eval(funcstrdict["u2"], xvectorlistsdictc)
        centerksu3=eval(funcstrdict["u3"], xvectorlistsdictc)
        print (r1, r2, r3)

        resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 1000, nvoly=1)
        mkprc=ex3.makepercent(leftksu2, rightksu2, leftksu3, rightksu3, resdict, draw=0)
        print ("Percent of good products: ", mkprc[0])


        if max(bestpercent, mkprc[0])!=bestpercent:
            bestpercent=max (bestpercent, mkprc[0])
            vectr=r1,r2,r3

    print (bestpercent, "\n", vectr)


    #отображаем лучший вариант
    xvectorlistsdict = {"u1":[100],  "r1":[r1], "r2":[r2], "r3":[r3]}
    resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 10000, nvoly=1)
    mkprc=ex3.makepercent(leftksu2, rightksu2, leftksu3, rightksu3, resdict, draw=1)



#test area:

funcstrdict= {"u2":"u1* (r2+r3)/(r1+r2+r3)", "u3":"u1* r3/(r1+r2+r3)"}


spreadvarslist  = ["r1", "r2", "r3"]
V=np.array       ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])



mntkrleq()


