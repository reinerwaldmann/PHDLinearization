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

        if (leftksu2 < ksu2 <rightksu2) and (leftksu3 < ksu3 < rightksu3):
            numOfOk+=1
            seqksu2OK.append(ksu2)
            seqksu3OK.append(ksu3)




    #рисование


    if (draw):

        plt.ion()
        fig = plt.figure()
        plt.xlabel('seqksu2')
        plt.ylabel('seqksu3')
      #  plt.axis([0.94, 0.97, 0.86, 0.92])
        plt.plot(seqksu2, seqksu3, 'ro') # Returns a tuple of line objects, thus the comma

        lplots=[plt.plot ([leftksu2, leftksu2], [leftksu3, rightksu3]),
        plt.plot ([rightksu2, rightksu2], [leftksu3, rightksu3]),
        plt.plot ([leftksu2, rightksu2], [leftksu3, leftksu3]),
        plt.plot ([leftksu2, rightksu2], [rightksu3, rightksu3])]

        for p in lplots:
            plt.setp(p, color='b', linewidth=2.0)
        plt.draw()
        plt.show(block=True)




    return 100*numOfOk/numTotal, np.corrcoef(seqksu2, seqksu3), seqksu2, seqksu3
#



#test area

def ex3test():
    funcstrdict= {"u2":"u1* (r2+r3)/(r1+r2+r3)", "u3":"u1* r3/(r1+r2+r3)"}
    xvectorlistsdict = {"u1":[100],  "r1":[20], "r2":[30], "r3":[400]}
    spreadvarslist  = ["r1", "r2", "r3"]
    V=np.array       ( [[4, 2, 3],
                                    [2, 9, 6],
                                    [3, 6, 16]])

    resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 10000, nvoly=1)

    ksu2=0.955
    ksu3=0.888

    diap=0.001
    leftksu2=ksu2-diap
    rightksu2=ksu2+diap
    leftksu3=ksu3-diap
    rightksu3=ksu3+diap
    print (makepercent(leftksu2, rightksu2, leftksu3, rightksu3, resdict, draw=1))








#first - y axiss, second - x axis


#xvectorlistsdictc = {"u1":100,  "r1":20, "r2":30, "r3":400}
#print ("center of the ellise\n",eval(funcstrdict["u2"], xvectorlistsdictc), "\n",  eval(funcstrdict["u3"], xvectorlistsdictc) )