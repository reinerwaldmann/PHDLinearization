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
        plt.figure(1) # Here's the part I need, but numbering starts at 1!
        plt.xlabel('seqksu2')
        plt.ylabel('seqksu3')
        #plt.axis([0.94, 0.97, 0.86, 0.92])
        plt.plot(seqksu2, seqksu3, 'ro')
        plt.draw()



#        time.sleep(1)



        plt.plot(seqksu2OK, seqksu3OK, 'ro')
        plt.draw()


        plt.show()

        # plt.figure(2) # Here's the part I need, but numbering starts at 1!
        # #plt.axis([0.94, 0.97, 0.86, 0.92])
        # plt.xlabel('seqksu2')
        # plt.ylabel('seqksu3')
        # plt.plot(seqksu2OK, seqksu3OK, 'ro')
        # plt.draw()





    #plt.axis([0.8,0.97 ,0.8,0.97])
    #plt.figure(2) # Here's the part I need, but numbering starts at 1!



    return 100*numOfOk/numTotal, np.corrcoef(seqksu2, seqksu3)


#

# plt.plot([1,2,3,4], [1,4,9,16], 'ro')
# plt.axis([0, 6, 0, 20])
# plt.show()






#print (makepercent(0.95, 0.96, 0.89, 0.9, resdict ))

#for i in range (0,10):



funcstrdict= {"u2":"u1* (r2+r3)/(r1+r2+r3)", "u3":"u1* r3/(r1+r2+r3)"}
xvectorlistsdict = {"u1":[100],  "r1":[20], "r2":[30], "r3":[400]}
spreadvarslist  = ["r1", "r2", "r3"]
V=np.array       ( [[4, 2, 3],
                                [2, 9, 6],
                                [3, 6, 16]])

resdict=sg.generate (funcstrdict, xvectorlistsdict, spreadvarslist, V, 10000, nvoly=1)
print (makepercent(0.88, 0.90, 0.95, 0.96, resdict, draw=1))

#first - y axiss, second - x axis


#xvectorlistsdictc = {"u1":100,  "r1":20, "r2":30, "r3":400}
#print ("center of the ellise\n",eval(funcstrdict["u2"], xvectorlistsdictc), "\n",  eval(funcstrdict["u3"], xvectorlistsdictc) )