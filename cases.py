__author__ = 'reiner'

import math

#import matplotlib.pyplot as plt
import numpy as np
import sympy

from scipy import optimize


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

#    print (jac)
#    exit(0)

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
    Vbe=VIN-VCC
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

#    print (jac)
#    exit(0)

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
    IS=b[2]

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

    Is=b[2]

    jac[0][1]=-Is*(math.exp(Vbc/FT) - 1)/b[1]**2 - GMIN*Vbc/b[1]**2
    jac[0][0]=-Is*(math.exp(Vbe/FT) - 1)/b[0]**2 - GMIN*Vbe/b[0]**2
    jac[0][2]=(math.exp(Vbc/FT) - 1)/b[1] + (math.exp(Vbe/FT) - 1)/b[0]
    jac[0][3]=0

    jac[1][0]=0
    jac[1][1]=Is*(math.exp(Vbc/FT) - 1)/b[1]**2 + GMIN*Vbc/b[1]**2
    jac[1][2]=(1 - Vbc/b[3])*(-math.exp(Vbc/FT) + math.exp(Vbe/FT)) - (math.exp(Vbc/FT) - 1)/b[1]
    jac[1][3]=Is*Vbc*(-math.exp(Vbc/FT) + math.exp(Vbe/FT))/b[3]**2 + GMIN*Vbc*(-Vbc + Vbe)/b[3]**2

    return jac


def testEstimateErlie():

    #за 14-16 итераций с единичного вектора к результирующему с нулевой ошибкой, несмотря на разные порядки!!! при дисперсии 0.1 на 50 точках при обычном плане (учитывать большее число точек на самом деле)
    #эпичная победа



    """
    Пробуем произвести экстракцию параметров модели по параметрам транзистора Эрли
    :return:
    """

    jacf=lambda x,b,c,y: outTransParamErlieFormatJAC (x,b)
    funcf=lambda x,b,c: outTransParamErlieFormat (x,b)

    c={}
    Ve=np.array([ [0.1, 0],
                     [0, 0.1] ]  )
    #BF,BR,IS,VA
    #коэфф передачи по току в схеме с оэ нормальный режим, -//- реверсный, ток утечки, напряжение Эрли в активном режиме
    btrue=[120,1,1.28e-15, 10]
    #binit=[115,0.1,1, 11]

    binit=[1]*len(btrue)

    bstart=[100,0.5,1, 5]
    bend=[125,2,2, 15]




    xstart=[0.001,0.001]
    xend=[4,4]

    N=50

    print("performing normal research:")
    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    #надо добавить скажем априорный план, с фильтрованием точек


    optimized_measdata=o_p.filterList(measdata, lim=1e55)
    print ('Plan optimization: measdatalen={0} optimized={1}'.format(len(measdata), len(optimized_measdata )))
    measdata = optimized_measdata


    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))

    print (gknu[0])



    #aprior plan
    print("Performing aprior plan:")
    oplan = o_ap.grandApriornPlanning(xstart, xend, 10, bstart, bend, c, Ve, jacf, Ntries=5)
    o_p.writePlanToFile(oplan, 'Aprior plan Erlie')
    measdata = o_p.makeMeasAccToPlan(funcf, oplan, btrue, c, Ve)

    optimized_measdata=o_p.filterList(measdata, lim=1e55)
    print ('Plan optimization: measdatalen={0} optimized={1}'.format(len(measdata), len(optimized_measdata )))
    measdata = optimized_measdata


    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))




def testSquareFunction():
#проходит за 9 итераций b=b+deltab*mu
    xstart=[1, 100]
    xend=[20,200]
    N=10
    c={"a":1000}
    funcf=lambda x,b,c=None:  [ b[0]+b[1]*x[0]+b[2]*x[1]+b[3]*x[0]*x[1]+b[4]*x[0]**2+b[5]*x[1]**2,   b[6]+b[7]*x[0]+b[8]*x[1]+b[9]*x[0]*x[1]+b[10]*x[0]**2+b[11]*x[1]**2 ]
    jacf = lambda x,b,c=None,y=None: np.array([ [1, x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1], 0, 0, 0, 0, 0, 0],
                                       [0,0,0,0,0,0,1,x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1]] ])
    #убрал .T
    Ve=np.array([ [0.0001, 0],
                  [0, 0.0001]]  )
    bstart=[0.8,0.4,1.4,0.2,0.9,0.3,1.4,0.2,2.1,3.1,4.1,5.1]
    btrue=  [1.1,0.6,1.6,0.4,1.1,0.6,1.6,0.4,2.5,3.3,4.6,5.6]

    bend=list(np.array(btrue)+np.array(btrue)-np.array(bstart))
    binit=[1]*len(btrue)

    N=50
    # print("performing normal research:")
    # startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    # measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    # gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c,NSIG=10)
    # print (gknu)
    # print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))
    # print (gknu[0])

    print("performing aprior research:")
    oplan = o_ap.grandApriornPlanning (xstart, xend, 10, bstart, bend, c, Ve, jacf, func=None, Ntries=1)[1]
    print(oplan)

    measdata = o_p.makeMeasAccToPlan(funcf, oplan, btrue, c, Ve)
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c,NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))
    print (gknu[0])



    # plan=makeUniformExpPlan(xstart, xend, N)
    # func = lambda x,b,c: [x[0]*b[0]+c["a"], x[1]*b[1]+c["a"], x[2]*b[2]+c["a"]]
    # meas = makeMeasAccToPlan(func, plan,  b, c, [0.0001]*3)
    # for x in meas:
    #     print (x)

def testSimpleFunction ():
    #проходит, 2 итерации, если  b=b+deltab*mu


    funcstrdict= {"y1":"u1* (r2+r3)", "y2":"u1* r3"}
    xstart=[1, 100]
    xend=[20,200]
    N=50
    funcf=lambda x,b,c=None:  [x[1]*x[0]*(b[0]+b[1]), x[0]*b[1]+x[1]*b[0]]
    jacf = lambda x,b,c=None,y=None: np.array([[x[0]*x[1],x[0]*x[1]],[x[1], x[0]]])
    Ve=np.array([ [0.00001, 0],
                  [0, 0.00001]]  )
    btrue=  [20,30]
    binit=[1]*len(btrue)
    print("performing normal research:")
    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)
    measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, None, Ve)
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, None,NSIG=10)
    print (gknu)
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, None))
    print (gknu[0])






def testModel():
    rng=np.arange(0.01,2,0.01)
    #print (outTransParamWErlie(1,2))

    #снимем входную ВАХ
    resrng=[outTransParam(x,5)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    #снимем выходную ВАХ выходит непохоже на PSPICE, но похоже на учебник
    #resrng=[outTransParam(0.2,x)[1] for x in rng] # изменяем напряжение на коллекторе при постоянном напряжении на базе - снимаем ток коллектора.

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
    I=Is*(math.exp(Vd/(FT*N)) -1)







#Для экстракции параметров
testEstimateErlie()

#Для отображения графика
#testModel()


