__author__ = 'vasilev_is'

#first of all, explicit scheme, no resistors

import math
import matplotlib.pyplot as plt
import numpy as np



FT=0.02586419 #подогнанное по pspice


import math as mpm

#mpm.mp.dps=40



def funcGP (x,b):
    """
    Модель Гуммеля-Пуна
    :param x:
    :param b:
    :return:
    """

    global FT

    #pr = lambda  x: mpm.mpf(str(x))
    pr=float

    b = list(map(pr  ,b))
    x = list(map(pr ,x))


    Vbe = -1*x[0]
    Vbc = -1*x[1]




    IS = b[0]       # сквозной ток насыщения   1e-16
    BF = b[1]       # максимальное значение нормального коэфф усиления по току с ОЭ 100
    VAF = b[2]      # прямое напряжение Эрли    inf
    VAR = b[3]      # инверсное напряжение Эрли inf
    IKF =  b[4]     # ток перехода к высококу уровню инжекции inf
    ISE = b[5]      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
    NE = b[6]       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
    NR = b[7]       # коэфф неидеальности для диффузного тока в инв режиме  1
    NF = b[8]       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
    NC = b[9]       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
    BR = b[10]      # инверсный коэфф усиления тока в схеме с ОЭ 1
    IKR = b[11]     # ток перехода к высокому уровню инжекции в инверсном включении inf
    ISC = b[12]     # генерационно-рекомбинационный ток насыщения колекторного перех 0
    RE = b[13]      # всякоразные сопротивления
    RC = b[14]
    RB = b[15]

    # главные составляющие
    Qfrac = .5*(1+Vbc/VAF+Vbe/VAR)+(.25*(1+Vbc/VAF+Vbe/VAR)**2+(IS/IKF)*(mpm.exp(Vbe/(NF*FT))-1)+(IS/IKR)*(mpm.exp(Vbc/(NR*FT))-1 )   )**.5

    Qfrac = 1/Qfrac

    Icc = Qfrac * IS * (mpm.exp(Vbe/(NF*FT))-mpm.exp(Vbc/(NR*FT) ))



    Ibe = (IS/BF) * (mpm.exp(Vbe/(NF*FT))-1)
    Ibc = (IS/BR) * (mpm.exp(Vbc/(NR*FT))-1)

    # генерационно-рекомбинационные составляющие
    Ire = ISE*(mpm.exp(Vbe/(NE*FT))-1)



    Irc = ISC*(mpm.exp(Vbc/(NC*FT))-1)

    if (x[1]<0):
        Irc=0
        Ibc*=2



    Ie = Icc+Ibe+Ire

    Ic = Ibc+Irc-Icc
    #print (Ibc,Irc,Icc)


    Ib = Ie-Ic
    y=[Ie, Ic, Ib]
    return y



def funcEM (x,b):
    """
    Модель Гуммеля-Пуна
    :param x:
    :param b:
    :return:
    """

    b = list(map(float,b))
    x = list(map(float,x))

    Vbe = -1*x[0]
    Vbc = -1*x[1]




    IS = b[0]       # сквозной ток насыщения   1e-16
    BF = b[1]       # максимальное значение нормального коэфф усиления по току с ОЭ 100
    VAF = b[2]      # прямое напряжение Эрли    inf
    VAR = b[3]      # инверсное напряжение Эрли inf
    IKF =  b[4]     # ток перехода к высококу уровню инжекции inf
    ISE = b[5]      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
    NE = b[6]       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
    NR = b[7]       # коэфф неидеальности для диффузного тока в инв режиме  1
    NF = b[8]       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
    NC = b[9]       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
    BR = b[10]      # инверсный коэфф усиления тока в схеме с ОЭ 1
    IKR = b[11]     # ток перехода к высокому уровню инжекции в инверсном включении inf
    ISC = b[12]     # генерационно-рекомбинационный ток насыщения колекторного перех 0
    RE = b[13]      # всякоразные сопротивления
    RC = b[14]
    RB = b[15]




    Icc = IS * (math.exp(Vbe/(NF*FT))-math.exp(Vbc/(NR*FT) ))


    Ibe = (IS/BF) * (math.exp(Vbe/(NF*FT))-1)
    Ibc = (IS/BR) * (math.exp(Vbc/(NR*FT))-1)
    # генерационно-рекомбинационные составляющие



    Ie = Icc+Ibe

    Ic = Ibc+Icc



    Ib = Ie-Ic
    y=[Ie, Ic, Ib]
    return y


def func (x,b):
    global FT

    Vbe = x[0]
    Vbc = x[1]



    IS = b[0]       # сквозной ток насыщения   1e-16
    BF = b[1]       # максимальное значение нормального коэфф усиления по току с ОЭ 100
    VAF = b[2]      # прямое напряжение Эрли    inf
    VAR = b[3]      # инверсное напряжение Эрли inf
    IKF =  b[4]     # ток перехода к высококу уровню инжекции inf
    ISE = b[5]      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
    NE = b[6]       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
    NR = b[7]       # коэфф неидеальности для диффузного тока в инв режиме  1
    NF = b[8]       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
    NC = b[9]       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
    BR = b[10]      # инверсный коэфф усиления тока в схеме с ОЭ 1
    IKR = b[11]     # ток перехода к высокому уровню инжекции в инверсном включении inf
    ISC = b[12]     # генерационно-рекомбинационный ток насыщения колекторного перех 0
    RE = b[13]      # всякоразные сопротивления
    RC = b[14]
    RB = b[15]




    Is=IS
    Va=VAF
    Is /= (1+(Vbc/Va))

    Ict = Is * (math.exp(Vbe/(FT*NE)) - math.exp(Vbc/(FT*NC)))
    Ie = (Is/BF) * (math.exp(Vbe/(FT*NE)) - 1) + Ict
    Ic =  -1*(Is/BR) * (math.exp(Vbc/(FT*NC)) - 1)+ Ict

    Ib = -1*(Ie + Ic)

    y=[-1*Ie, Ic, Ib]

    # вот насчёт коэффициентов неидеальности вопросы, к сожалению


    return y





def test ():
# .MODEL KT315 NPN (IS=10F BF=584.517 VAF=100 IKF=29.2714M ISE=131.803P
# + NE=2.08337 BR=1.95214 IKR=9.99996M ISC=100.316P RE=1 RC=5.48635
# + CJE=27.3893P VJE=700.001M MJE=500.287M CJC=27.3893P VJC=700.001M
# + MJC=500.287M TF=450.287P XTF=499.984M VTF=10 ITF=10.2268M TR=153.383P)

    inf=10e10

    IS = 10e-15
    BF = 584.517
    VAF = 100
    VAR = inf
    IKF =  29.2714e-3     # ток перехода к высококу уровню инжекции inf
    ISE = 131.803e-12      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
    NE = 2.08337       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
    NR = 1       # коэфф неидеальности для диффузного тока в инв режиме  1
    NF = 1       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
    NC = 1       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
    BR = 1.95214      # инверсный коэфф усиления тока в схеме с ОЭ 1
    IKR = 9.99996e-3     # ток перехода к высокому уровню инжекции в инверсном включении inf
    ISC = 100.316e-12     # генерационно-рекомбинационный ток насыщения колекторного перех 0
    RE = 1      # сопротивления эмиттера, коллектора, базы 0 0 0
    RC = 5.48635
    RB = 0



    b = [IS,
        BF,
        VAF,
        VAR,
        IKF,
        ISE,
        NE,
        NR,
        NF,
        NC,
        BR,
        IKR,
        ISC,
        RE,
        RC,
        RB,
]

    #1 input characteristics, common base, input on emitter

    Vbc=0

    v1 = -2.0
    v2 = 2.0

    fxed=.9


    funcw = funcGP

    # ВНИМАНИЕ! Модель работает только в нормальном активном режиме и в режиме отсечки
    # то есть только при положительных напряжениях, поданных на коллектор. При отрицательных кажет бред.



    print ([fxed,fxed], funcw([fxed,fxed],b))
#    print ([-fxed,-fxed], funcw([-fxed,-fxed],b))
 #   print ([fxed,-fxed], funcw([fxed,-fxed],b))
    print ([-fxed,fxed], funcw([-fxed,fxed],b))



    exit(0)


    xrange = np.arange (v1, v2, 0.00001)
    yrange1 = [funcw([Vbe, fxed] ,b)[0] for Vbe in xrange ]

    #plt.plot(xrange, yrange,  'r')
    plt.plot(xrange, yrange1,  'b')
    #plt.plot(xrange, yrange2,  'g')


    plt.ylabel('Ie')
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





    # попытка построить выходную закончилась фейлом, ибо там надо фиксировать ток
    # эмиттера, что пока непонятно, как сделать вообще

    # 2 output characteristics, common base, input on emitter


    yrange1 = [funcw([fxed, Vbc] ,b)[1] for Vbc in xrange ]


    plt.plot(xrange, yrange1,  'b')
    #plt.plot(xrange, yrange2,  'g')

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










