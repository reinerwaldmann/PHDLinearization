__author__ = 'vasilev_is'

#first of all, explicit scheme, no resistors

import math
import matplotlib.pyplot as plt
import numpy as np



FT=0.02586419 #подогнанное по pspice



def func (x,b):
    global FT

    Vbe = x[0]
    Vbc = x[1]

    Is = b[0]
    BF = b[1]
    BR = b[2]

    N = b[3]
    Va = b[4]

    Is /= (1+(Vbc/Va))

    Ict = Is * ( math.exp(Vbe/(FT*N)) - math.exp(Vbc/(FT*N)))
    Ie = (Is/BF) * (math.exp(Vbe/(FT*N)) - 1) + Ict
    Ic =  -1*(Is/BR) * (math.exp(Vbc/(FT*N)) - 1)+ Ict

    Ib = -1*(Ie + Ic)

    y=[Ie, Ic, Ib]

    # вот насчёт коэффициентов неидеальности вопросы, к сожалению


    return y

def test ():
# .MODEL KT315 NPN (IS=10F BF=584.517 VAF=100 IKF=29.2714M ISE=131.803P
# + NE=2.08337 BR=1.95214 IKR=9.99996M ISC=100.316P RE=1 RC=5.48635
# + CJE=27.3893P VJE=700.001M MJE=500.287M CJC=27.3893P VJC=700.001M
# + MJC=500.287M TF=450.287P XTF=499.984M VTF=10 ITF=10.2268M TR=153.383P)

    b = [10.0e-15, 584.517, 1.95214, 1, 100]

    #1 input characteristics, common base, input on emitter

    Vbc=0
    xrange = np.arange (.001, 1.5, 0.01)
    yrange = [func([Vbe, .9] ,b)[0] for Vbe in xrange ]
    yrange1 = [func([Vbe, -.0000001] ,b)[0] for Vbe in xrange ]

    plt.plot(xrange, yrange,  'r')
    plt.plot(xrange, yrange1,  'b')
    #plt.plot(xrange, yrange2,  'g')


    plt.ylabel('Ic')
    plt.xlabel('Ube')
    plt.grid()
    plt.show()

    #
    # по анализу:
    # полученная ВАХ соответствует учебнику. Реально узкий пучок, и всё такое
    # важные моменты: 1. Уравнения, по крайней мере Ict и Ie - ВЕРНЫ, и проверяются по уравнениям Кирхгофа
    # в учебничке или где-то там ещё неправильно указано направление токов, это, разумеется, переворачивает ВАХ
    #
    # 2. Пучок в данном варианте реально узкий. Потому что напряжение Эрли не введено




    exit(0)

    # попытка построить выходную закончилась фейлом, ибо там надо фиксировать ток
    # эмиттера, что пока непонятно, как сделать вообще

    # 2 output characteristics, common base, input on emitter

    xrange = np.arange (.001, 5, 0.01)

    yrange1 = [func([2, Vbc] ,b)[1] for Vbc in xrange ]
    yrange2 = [func([-.2, Vbc] ,b)[1] for Vbc in xrange ]


    plt.plot(xrange, yrange1,  'b')
    plt.plot(xrange, yrange2,  'g')

    plt.ylabel('Ic')
    plt.xlabel('Ubc')
    plt.grid()
    plt.show()

    # как же нам построить ВАХ, если для этого надо зафиксировать
    # ток??
    # Проблема в том, что у нас напряжения - входные параметры, токи - выходные (шо кстати, можно б изменить,
    # задав ток нагрузкой)
    # А транзистор - токовый (!) прибор. Вот и как получить  зависимость Iс от Ucb при ПОСТОЯННОМ Ie
    # есть дурацкий способ отбора точек ))))



test()










