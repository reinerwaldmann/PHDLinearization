__author__ = 'reiner'

import math

#import matplotlib.pyplot as plt
import numpy
import sympy


def testTransistorEMModel ():
    #схема с общим эмиттером
    VIN=0.2
    VCC=5

    #входные напряжения
    Vbe=VIN
    Vbc=VIN-VCC

    #параметры (пока не в стандарте)
    BF=120 #коэфф передачи по току в схеме с оэ нормальный режим
    BR=1 #коэфф передачи по току в схеме с оэ инверсный режим
    IS=1.28e-15 #ток утечки
    FT=0.026 #тепловой потенциал

    Ict=IS*(math.exp(Vbe/FT)-math.exp(Vbc/FT))
    Icc_sep_BF=(IS/BF)*(math.exp(Vbe/FT)-1)
    Iec_sep_BR=(IS/BR)*(math.exp(Vbc/FT)-1)

    Ic=Ict-Iec_sep_BR

    Ie=-Icc_sep_BF-Ict

    Ib=Icc_sep_BF+Iec_sep_BR

    print (Ic, Ie, Ib)

#testTransistorEMModel()


def outTransParam (Vin, Vcc):
    #схема с общим эмиттером
    VIN=Vin
    VCC=Vcc

    #входные напряжения
    Vbe=VIN
    Vbc=VIN-VCC

    #параметры (пока не в стандарте)
    BF=120 #коэфф передачи по току в схеме с оэ нормальный режим
    BR=1 #коэфф передачи по току в схеме с оэ инверсный режим
    IS=1.28e-15 #ток утечки
    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm


    Ict=IS*(math.exp(Vbe/FT)-math.exp(Vbc/FT))

    Icc_sep_BF=(IS/BF)*(math.exp(Vbe/FT)-1)

    Iec_sep_BR=(IS/BR)*(math.exp(Vbc/FT)-1)

    Ic=Ict-Iec_sep_BR

    Ie=-Icc_sep_BF-Ict

    Ib=Icc_sep_BF+Iec_sep_BR

    return  Ib, Ic, Ie



def outTransParamWErlie (Vin, Vcc):
    #Схема с общим эмиттером
    #Работает, выдаёт ВАХ, похожие на правду. Теперь надо, чтоб хоть какая-то из моих моделей полностью билась с PSPICE

    VIN=Vin
    VCC=Vcc

    #входные напряжения
    Vbe=VIN
    Vbc=VIN-VCC

    #параметры (пока не в стандарте)
    BF=120 #коэфф передачи по току в схеме с оэ нормальный режим
    BR=1 #коэфф передачи по току в схеме с оэ инверсный режим
    IS=1.28e-15 #ток утечки

    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm

    ISZERO=IS
    VA=1e1 #напряжение Эрли по умолчанию бесконечность

    GMIN=1e-12
    Ict=(ISZERO/(1+Vbc/VA)) *(math.exp(Vbe/FT)-math.exp(Vbc/FT)) #VA=INF даёт нам обычное выражение

    Ic=IS*( (math.exp(Vbe/FT)-math.exp(Vbc/FT))*(1-Vbc/VA)-(1/BR)*(math.exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/VA)-Vbc/BR)
    Ib=IS*( (1/BF)*(math.exp(Vbe/FT)-1)+(1/BR)*(math.exp(Vbc/FT)-1)) + GMIN*(Vbe/BF+Vbc/BR)

    return  Ib, Ic

def outTransParamErlieFormat (x,b,c):
    """
    :param x: вектор входов (напряженения  база, коллектор)
    :param b: вектор коэффициентов: BF,BR,IS,VA
    :param c: пустой, фиксированные параметры заданы в функции
    :return: ток базы, ток коллектора
    """

    VIN=x[0]
    VCC=x[1]

    #входные напряжения
    Vbe=VIN
    Vbc=VIN-VCC

    #параметры (пока не в стандарте)
    BF=120 #коэфф передачи по току в схеме с оэ нормальный режим
    BR=1 #коэфф передачи по току в схеме с оэ инверсный режим
    IS=1.28e-15 #ток утечки

    BF=b[0]
    BR=b[1]
    IS=b[2]

    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm
    ISZERO=IS

    VA=1e1 #напряжение Эрли по умолчанию бесконечность
    VA=b[3]

    GMIN=1e-12

#    Ict=(ISZERO/(1+Vbc/VA)) *(math.exp(Vbe/FT)-math.exp(Vbc/FT)) #VA=INF даёт нам обычное выражение
    Ic=IS*( (math.exp(Vbe/FT)-math.exp(Vbc/FT))*(1-Vbc/VA)-(1/BR)*(math.exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/VA)-Vbc/BR)
    Ib=IS*( (1/BF)*(math.exp(Vbe/FT)-1)+(1/BR)*(math.exp(Vbc/FT)-1)) + GMIN*(Vbe/BF+Vbc/BR)

    return  Ib, Ic



def outTransParamErlieFormatJAC (x,b,c):
    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm
    GMIN=1e-12

     #входные напряжения
    VIN=x[0]
    VCC=x[1]
    Vbe=VIN
    Vbc=VIN-VCC

    jac=numpy.zeros((2, 3))

    jac[0][0]=-b[2]*(math.exp(Vbe/FT) - 1)/b[0]**2 - GMIN*Vbe/b[0]**2
    jac[0][1]=-b[2]*(math.exp(Vbc/FT) - 1)/b[1]**2 - GMIN*Vbc/b[1]**2
    jac[0][2]=(math.exp(Vbc/FT) - 1)/b[1] + (math.exp(Vbe/FT) - 1)/b[0]

    jac[1][0]=0
    jac[1][1]=b[2]*(math.exp(Vbc/FT) - 1)/b[1]**2 + GMIN*Vbc/b[1]**2
    jac[1][2]=(1 - Vbc/b[3])*(-math.exp(Vbc/FT) + math.exp(Vbe/FT)) - (math.exp(Vbc/FT) - 1)/b[1]

    return jac



def getderivative_outTransParamErlieFormatJAC():
    """
    Якобиан по коэфффициентам
    """
    funstr=["B2*( (1/B0)*(exp(Vbe/FT)-1)+(1/B1)*(exp(Vbc/FT)-1)) + GMIN*(Vbe/B0+Vbc/B1)", "B2*( (exp(Vbe/FT)-exp(Vbc/FT))*(1-Vbc/B3)-(1/B1)*(exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/B3)-Vbc/B1)" ]

    resstr=""

    for i in range (0, len(funstr)):
        for ind in range (0, 3):
            #print(sympy.diff(funstr[i],'B{0}'.format(ind)))


            resstr+=sympy.diff(funstr[i],'B{0}'.format(ind)).__str__()
            resstr+="\n"
        resstr+="------------------\n"

    return resstr



print(outTransParamErlieFormatJAC (None,None,None))
exit(0)




def outTransParamPopov (Vin, Vcc):
    #Схема с общей базой!!!!!!!!!


    W1=1 # напряжение на эмиттере
    W2=2 # напряжение на коллекторе

    VIN=Vin
    VCC=Vcc

    #входные напряжения
    Vbe=VIN
    Vbc=VIN-VCC

    #параметры (пока не в стандарте)
    BF=120 #коэфф передачи по току в схеме с оэ нормальный режим
    BR=1 #коэфф передачи по току в схеме с оэ инверсный режим
    ISC=1.55e-15 #ток утечки
    FT=0.026 #тепловой потенциал

    ISZERO=IS
    VA=1e1 #напряжение Эрли по умолчанию бесконечность

    GMIN=1e-12
    Ict=(ISZERO/(1+Vbc/VA)) *(math.exp(Vbe/FT)-math.exp(Vbc/FT)) #VA=INF даёт нам обычное выражение

    Ic=IS*( (math.exp(Vbe/FT)-math.exp(Vbc/FT))*(1-Vbc/VA)-(1/BR)*(math.exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/VA)-Vbc/BR)
    Ib=IS*( (1/BF)*(math.exp(Vbe/FT)-1)+(1/BR)*(math.exp(Vbc/FT)-1)) + GMIN*(Vbe/BF+Vbc/BR)

    return  Ib, Ic




def testModel():
    rng=numpy.arange(0.01,2,0.01)
    #print (outTransParamWErlie(1,2))
    #TODO ну и три возможные задачи: 1. Можно дальше улучшать модель, доведя до Гуммеля-Пуна
    #TODO 2. Обязательно попробовать оценку через Гаусса-Ньютона
    #TODO 3. Сделать модель IGBT или ещё чего-нибудь такого
    #resrng=[outTransParam(x,0)[1] for x in rng] # изменяем напряжение на базе при постоянном напряжении на колллекторе - снимаем ток коллектора.
    #в каком
    #снимем входную ВАХ
    resrng=[outTransParam(x,5)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на колллекторе - снимаем ток коллектора.
    #снимем выходную ВАХ выходит непохоже на PSPICE, но похоже на учебник
    #resrng=[outTransParam(0.2,x)[1] for x in rng] # изменяем напряжение на базе при постоянном напряжении на колллекторе - снимаем ток коллектора.

    #resrng=[outTransParam(0,x)[1] for x in rng]    #[1] - это ток коллектора
    # resrng1=[outTransParam(2,x)[1] for x in rng]
    #resrng=[outTransParamWErlie(0,x)[1] for x in rng]
    #resrng1=[outTransParamWErlie(2,x)[1] for x in rng]
    plt.plot(rng , resrng)
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()
    #ax.set_xticks(numpy.arange(0,1,0.1))
    #ax.set_yticks(numpy.arange(0,1.,0.1))
    #print (resrng1)
    #print (resrng)
    #plt.plot(rng , resrng1)
    plt.show()


