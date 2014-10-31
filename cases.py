__author__ = 'reiner'

import math

#import matplotlib.pyplot as plt
import numpy as np
import sympy

import Ofiura_Estimation as o_e
import Ofiura_planning as o_p
import Ofiura_Qualitat as o_q
import Ofiura_ApriorPlanning as o_ap



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

def outTransParamFormat (x,b,c=None):
     #схема с общим эмиттером
    VIN=x[0]
    VCC=x[1]

    #входные напряжения
    Vbe=VIN
    Vbc=VIN-VCC

    #параметры (пока не в стандарте)
    BF=120 #коэфф передачи по току в схеме с оэ нормальный режим
    BR=1 #коэфф передачи по току в схеме с оэ инверсный режим
    IS=1.28e-15 #ток утечки
    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm

    BF=b[0]
    BR=b[1]
    IS=b[2]*1e-15
    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm

    Ict=IS*(math.exp(Vbe/FT)-math.exp(Vbc/FT))

    Icc_sep_BF=(IS/BF)*(math.exp(Vbe/FT)-1)
    Iec_sep_BR=(IS/BR)*(math.exp(Vbc/FT)-1)

    Ie=-Icc_sep_BF-Ict

    Ib=Icc_sep_BF+Iec_sep_BR
    Ic=Ict-Iec_sep_BR


    return  [Ib, Ic]

def getderivative_outTransParamFormatJAC():
    """
    Якобиан по коэфффициентам
    """

    funstr = ["(b2*1e-15/b0)*(exp(Vbe/FT)-1) + (b2*1e-15/b1)*(exp(Vbc/FT)-1)", "b2*1e-15*(exp(Vbe/FT)-exp(Vbc/FT)) - (b2*1e-15/b1)*(exp(Vbc/FT)-1)"]

    resstr=""

    for i in range (0, len(funstr)):
        for ind in range (0, 4):
            #print(sympy.diff(funstr[i],'B{0}'.format(ind)))


            resstr+=sympy.diff(funstr[i],'b{0}'.format(ind)).__str__()
            resstr+="\n"
            print ('b{0}'.format(ind))

        resstr+="------------------\n"

    return resstr



def outTransParamFormatJAC (x,b,c=None):

    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm
    GMIN=1e-12

     #входные напряжения
    VIN=x[0]
    VCC=x[1]
    Vbe=VIN
    Vbc=VIN-VCC

    jac=np.zeros((2, 3))


    #коэфф передачи по току в схеме с оэ нормальный режим, -//- реверсный, ток утечки


    jac[0][0]=-1.0e-15*b[2]*(math.exp(Vbe/FT) - 1)/b[0]**2
    jac[0][1]=-1.0e-15*b[2]*(math.exp(Vbc/FT) - 1)/b[1]**2
    jac[0][2]=1.0e-15*(math.exp(Vbc/FT) - 1)/b[1] + 1.0e-15*(math.exp(Vbe/FT) - 1)/b[0]

    jac[1][0]=0
    jac[1][1]=1.0e-15*b[2]*(math.exp(Vbc/FT) - 1)/b[1]**2
    jac[1][2]=-1.0e-15*math.exp(Vbc/FT) + 1.0e-15*math.exp(Vbe/FT) - 1.0e-15*(math.exp(Vbc/FT) - 1)/b[1]

    print (jac)
    exit(0)

    return jac

def testEstimate():
    """
    Пробуем произвести экстракцию параметров модели по параметрам транзистора Эрли
    :return:
    """
    jacf=lambda x,b,c,y: outTransParamFormatJAC (x,b)
    funcf=lambda x,b,c: outTransParamFormat (x,b)

    c={}
    Ve=np.array([ [0.000001, 0],
                     [0, 0.000001] ]  )

    #BF,BR,IS
    #коэфф передачи по току в схеме с оэ нормальный режим, -//- реверсный, ток утечки
    btrue=[120,1,1.28]
    binit=[110,2,1.28]

    bstart=[100,0.5,1]
    bend=[125,2,2]

    xstart=[0.001,0.001]
    xend=[5,5]

    N=50

    print("performing normal research:")
    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    #надо добавить скажем априорный план, с фильтрованием точек

    print ('Plan optimization: measdatalen={0} optimized={1}'.format(len(measdata), len(list(filter(for_filter, measdata)) )))
    measdata = list(filter(for_filter, measdata))


    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))

    print (gknu[0])

    return

#     aprior plan
    print("Performing aprior plan:")
    oplan = o_ap.grandApriornPlanning(xstart, xend, 10, bstart, bend, c, Ve, jacf, Ntries=2)
    o_p.writePlanToFile(oplan, 'Aprior_plan')
    measdata = o_p.makeMeasAccToPlan(funcf, oplan, btrue, c, Ve)
    filteredmeasdata=list(filter(for_filter, measdata))
    print ('Plan optimization: measdatalen={0} optimized={1}'.format(len(measdata), len(filteredmeasdata) ))



    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))







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

def outTransParamErlieFormat (x,b,c=None):
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
    IS=b[2]*1e-15

    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm
    ISZERO=IS

    VA=1e1 #напряжение Эрли по умолчанию бесконечность
    VA=b[3]

    GMIN=1e-12

#    Ict=(ISZERO/(1+Vbc/VA)) *(math.exp(Vbe/FT)-math.exp(Vbc/FT)) #VA=INF даёт нам обычное выражение

    Ib=IS*( (1/BF)*(math.exp(Vbe/FT)-1)+(1/BR)*(math.exp(Vbc/FT)-1)) + GMIN*(Vbe/BF+Vbc/BR)
    Ic=IS*( (math.exp(Vbe/FT)-math.exp(Vbc/FT))*(1-Vbc/VA)-(1/BR)*(math.exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/VA)-Vbc/BR)

    return  [Ib, Ic]

def outTransParamErlieFormatJAC (x,b,c=None):
    FT=0.026 #тепловой потенциал - он фиксирован для 27С http://cads.narod.ru/kurs/OrCAD.htm
    GMIN=1e-12

     #входные напряжения
    VIN=x[0]
    VCC=x[1]
    Vbe=VIN
    Vbc=VIN-VCC

    jac=np.zeros((2, 4))

    jac[0][0]=-b[2]*1e-15*(math.exp(Vbe/FT) - 1)/b[0]**2 - GMIN*Vbe/b[0]**2
    jac[0][1]=-b[2]*1e-15*(math.exp(Vbc/FT) - 1)/b[1]**2 - GMIN*Vbc/b[1]**2
    jac[0][2]=(math.exp(Vbc/FT) - 1)/b[1] + (math.exp(Vbe/FT) - 1)/b[0]
    jac[0][3]=0

    jac[1][0]=0
    jac[1][1]=b[2]*1e-15*(math.exp(Vbc/FT) - 1)/b[1]**2 + GMIN*Vbc/b[1]**2
    jac[1][2]=(1 - Vbc/b[3])*(-math.exp(Vbc/FT) + math.exp(Vbe/FT)) - (math.exp(Vbc/FT) - 1)/b[1]
    jac[1][3]=b[2]*1e-15*Vbc*(-math.exp(Vbc/FT) + math.exp(Vbe/FT))/b[3]**2 + GMIN*Vbc*(-Vbc + Vbe)/b[3]**2

    return jac

def getderivative_outTransParamErlieFormatJAC():
    """
    Якобиан по коэфффициентам
    """
    funstr=["B2*( (1/B0)*(exp(Vbe/FT)-1)+(1/B1)*(exp(Vbc/FT)-1)) + GMIN*(Vbe/B0+Vbc/B1)", "B2*( (exp(Vbe/FT)-exp(Vbc/FT))*(1-Vbc/B3)-(1/B1)*(exp(Vbc/FT)-1))+GMIN*((Vbe-Vbc)*(1-Vbc/B3)-Vbc/B1)" ]

    resstr=""

    for i in range (0, len(funstr)):
        for ind in range (0, 4):
            #print(sympy.diff(funstr[i],'B{0}'.format(ind)))


            resstr+=sympy.diff(funstr[i],'B{0}'.format(ind)).__str__()
            resstr+="\n"
            print ('B{0}'.format(ind))

        resstr+="------------------\n"

    return resstr



def for_filter (x):
    for val in x['y']:
        if val>1e55:
            return False
            print ("ex")
    return True


def testEstimateErlie():
    """
    Пробуем произвести экстракцию параметров модели по параметрам транзистора Эрли
    :return:
    """

    jacf=lambda x,b,c,y: outTransParamErlieFormatJAC (x,b)
    funcf=lambda x,b,c: outTransParamErlieFormat (x,b)

    c={}
    Ve=np.array([ [0.00001, 0],
                     [0, 0.00001] ]  )


    #BF,BR,IS,VA
    #коэфф передачи по току в схеме с оэ нормальный режим, -//- реверсный, ток утечки, напряжение Эрли в активном режиме
    btrue=[120,1,1.28, 10]
    binit=[1,1,1, 1]


    bstart=[100,0.5,1, 5]
    bend=[125,2,2, 15]



#РАЗНЫЙ ПОРЯДОК КОЭФФИЦИЕНТОВ, вот что может всё портить!!!!

    xstart=[0.001,0.001]
    xend=[4,4]

    N=50

    print("performing normal research:")
    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    #надо добавить скажем априорный план, с фильтрованием точек

    print ('Plan optimization: measdatalen={0} optimized={1}'.format(len(measdata), len(list(filter(for_filter, measdata)) )))
    measdata = list(filter(for_filter, measdata))


    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))

    print (gknu[0])



    #aprior plan
    # print("Performing aprior plan:")
    # oplan = o_ap.grandApriornPlanning(xstart, xend, 10, bstart, bend, c, Ve, jacf, Ntries=5)
    # o_p.writePlanToFile(oplan, 'Aprior plan Erlie')
    # measdata = o_p.makeMeasAccToPlan(funcf, oplan, btrue, c, Ve)
    # filteredmeasdata=list(filter(for_filter, measdata))
    # print ('Plan optimization: measdatalen={0} optimized={1}'.format(len(measdata), len(filteredmeasdata) ))
    #
    #
    #
    # gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=10)
    # print (gknu)
    # print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))



testEstimate()




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


#Метод считает начальное приближение уже весьма хорошим (предлагает слишком маленький шаг)
#возможная причина - производные по коэффициентам слишком малы. (порядка 10^-19) относительно самих коэффициентов
#то есть выходные токи сильно больше зависят от напряжений на базе и коллекторе, нежели от каких-то там параметров.
#выходит, хорошо, когда выходные переменые примерно одинаково зависят от оцениваемых параметров и входных переменых.
