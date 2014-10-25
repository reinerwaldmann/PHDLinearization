__author__ = 'reiner'

import math

import matplotlib.pyplot as plt
import numpy


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
    FT=0.026 #тепловой потенциал

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
    FT=0.026 #тепловой потенциал

    ISZERO=IS
    VA=1e1 #напряжение Эрли по умолчанию бесконечность

    GMIN=1e-12
    Ict=(ISZERO/(1+Vbc/VA)) *(math.exp(Vbe/FT)-math.exp(Vbc/FT)) #VA=INF даёт нам обычное выражение

    Ic=IS*( (math.exp(Vbe/FT)-math.exp(Vbc/FT))*(1-Vbc/VA)-(1/BR)*(math.exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/VA)-Vbc/BR)
    Ib=IS*( (1/BF)*(math.exp(Vbe/FT)-1)+(1/BR)*(math.exp(Vbc/FT)-1)) + GMIN*(Vbe/BF+Vbc/BR)

    return  Ib, Ic



rng=numpy.arange(0.01,1,0.01)


#TODO ну и три возможные задачи: 1. Можно дальше улучшать модель, доведя до Гуммеля-Пуна
#TODO 2. Обязательно попробовать оценку через Гаусса-Ньютона
#TODO 3. Сделать модель IGBT или ещё чего-нибудь такого


# resrng=[outTransParam(x,50)[0] for x in rng]


# resrng=[outTransParam(0,x)[1] for x in rng]
# resrng1=[outTransParam(2,x)[1] for x in rng]

resrng=[outTransParamWErlie(0,x)[1] for x in rng]
resrng1=[outTransParamWErlie(2,x)[1] for x in rng]





plt.plot(rng , resrng1)
#plt.axis([0.0,1.0,-2,4])

plt.grid()
plt.show()


#ax.set_xticks(numpy.arange(0,1,0.1))
#ax.set_yticks(numpy.arange(0,1.,0.1))




print (resrng1)



#print (resrng)
#plt.plot(rng , resrng1)




plt.show()


