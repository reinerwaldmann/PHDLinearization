__author__ = 'reiner'
import matplotlib.pyplot as plt
import numpy as np

import sequence_generation as sg


#Этот файл и функции, в нём определённые, относятся только к частной задаче исследования делителя напряжения о трёх резисторах.


def makepercent (leftksu2, rightksu2, leftksu3, rightksu3, generatedvals, draw=1):
    """
    Выполняет задачу 3 исследования делителей. Принимает на вход сгенерированные значения, включая выходные и входные напряжения.
    Делит выходные на входные, таким образом получая коэффициенты деления. Сличает их с диапазонами rangeu2 и rangeu3, на основе этого
    принимает решение, годен делитель или же не годен. Если годен, прибавляет число годных. если нет - прибавляет число негодных.

    Рисует точки, отображающие коэфф. деления. Они должны составить собою эллипс.
    """
    numTotal=len(generatedvals)
    numOfOk=0

    seqksu2=list()
    seqksu3=list()

    seqksu2OK=list()
    seqksu3OK=list()

    for set in generatedvals:
        ksu2=set["u2"]/set["u1"]
        ksu3=set["u3"]/set["u1"]

        seqksu2.append(ksu2)
        seqksu3.append(ksu3)

        if (leftksu2 < ksu3 <rightksu2) and (leftksu3 < ksu2 < rightksu3):
            numOfOk+=1
            seqksu2OK.append(ksu2)
            seqksu3OK.append(ksu3)




    #рисование


    if (draw):
        pass
        # plt.ion()
        # fig = plt.figure()
        # plt.xlabel('seqksu2')
        # plt.ylabel('seqksu3')
        # plt.axis([0.94, 0.97, 0.86, 0.92])
        # plt.plot(seqksu2, seqksu3, 'ro') # Returns a tuple of line objects, thus the comma
        #
        # lplots=[plt.plot ([leftksu3, leftksu3], [leftksu2, rightksu2]),
        # plt.plot ([rightksu3, rightksu3], [leftksu2, rightksu2]),
        # plt.plot ([leftksu3, rightksu3], [leftksu2, leftksu2]),
        # plt.plot ([leftksu3, rightksu3], [rightksu2, rightksu2])]
        #
        # for p in lplots:
        #     plt.setp(p, color='b', linewidth=2.0)




        # plt.draw()
        # time.sleep(10)
        #
        # plt.clf()
        #
        #
        # plt.axis([0.94, 0.97, 0.86, 0.92])
        # plt.plot(seqksu2OK, seqksu3OK, 'ro') # Returns a tuple of line objects, thus the comma
        #
        #
        # plt.draw()
        # time.sleep(1)
        #plt.show(block=True)




    return 100*numOfOk/numTotal, np.corrcoef(seqksu2, seqksu3), seqksu2, seqksu3






#Monte-Karlo method V LOB:

bestpercent=0

vectr=[]

for r1 in range(10,500, 1):
    for r2 in range(10,500, 1):
        for r3 in range(10,500, 1):
            funcstrdict= {"u2":"u1* (r2+r3)/(r1+r2+r3)", "u3":"u1* r3/(r1+r2+r3)"}
            xvectorlistsdict = {"u1":[100],  "r1":[r1], "r2":[r2], "r3":[r3]}
            xvectorlistsdictc = {"u1":100,  "r1":r1, "r2":r2, "r3":r3}


            spreadvarslist  = ["r1", "r2", "r3"]
            V=np.array       ( [[4, 2, 3],
                                    [2, 9, 6],
                                    [3, 6, 16]])

            leftksu2, rightksu2, leftksu3, rightksu3 = 0.88, 0.90, 0.95, 0.96

            centerksu2=eval(funcstrdict["u2"], xvectorlistsdictc)
            centerksu3=eval(funcstrdict["u3"], xvectorlistsdictc)


            centerksu2req=leftksu2+(rightksu2-leftksu2)
            centerksu3req=leftksu3+(rightksu3-leftksu3)

            border=10


            if centerksu2req-border<centerksu2<centerksu2req+border and centerksu3req-border<centerksu3<centerksu3req+border:
                pass

            else:
                continue

            print (centerksu2, centerksu3)

            resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 1000, nvoly=1)

            mkprc=makepercent(leftksu2, rightksu2, leftksu3, rightksu3, resdict, draw=0)
            print ("Percent of good products: ", mkprc[0])


            if max(bestpercent, mkprc[0])!=bestpercent:
                bestpercent=max (bestpercent, mkprc[0])
                vectr=r1,r2,r3



print ("bestpercent=",bestpercent)
print ("vectr=",vectr)



exit(0)










plt.ion()
fig = plt.figure()

if 0:
    funcstrdict= {"u2":"u1* (r2+r3)/(r1+r2+r3)", "u3":"u1* r3/(r1+r2+r3)"}
    #xvectorlistsdict = {"u1":[100],  "r1":[20], "r2":[r2], "r3":[400]}
    xvectorlistsdict = {"u1":[100],  "r1":[1538], "r2":[920], "r3":[12193]}
    spreadvarslist  = ["r1", "r2", "r3"]
    V=np.array       ( [[4, 2, 3],
                                    [2, 9, 6],
                                    [3, 6, 16]])

    resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 10000, nvoly=1)
    leftksu2, rightksu2, leftksu3, rightksu3 = 0.88, 0.90, 0.95, 0.96
    mkprc=makepercent(leftksu2, rightksu2, leftksu3, rightksu3, resdict, draw=1)

    #print ("r2=",r2)
    print ("Percent of good products: ", mkprc[0])
    print ("CorrCoef: ", mkprc[1][0][1])

#    if mkprc[1][0][1]<0:
 #       break



    plt.xlabel('seqksu2')
    plt.ylabel('seqksu3')
    #plt.axis([0.94, 0.97, 0.86, 0.92])


    plt.plot(mkprc[2], mkprc[3], 'ro') # Returns a tuple of line objects, thus the comma

    lplots=[plt.plot ([leftksu3, leftksu3], [leftksu2, rightksu2]),
    plt.plot ([rightksu3, rightksu3], [leftksu2, rightksu2]),
    plt.plot ([leftksu3, rightksu3], [leftksu2, leftksu2]),
    plt.plot ([leftksu3, rightksu3], [rightksu2, rightksu2])]

    for p in lplots:
        plt.setp(p, color='b', linewidth=2.0)

    plt.draw()
    #time.sleep(1)

#    plt.clf()

plt.show(block=True)





#first - y axiss, second - x axis


#xvectorlistsdictc = {"u1":100,  "r1":20, "r2":30, "r3":400}
#print ("center of the ellise\n",eval(funcstrdict["u2"], xvectorlistsdictc), "\n",  eval(funcstrdict["u3"], xvectorlistsdictc) )