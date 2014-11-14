__author__ = 'reiner'
import math

#import matplotlib.pyplot as plt
import numpy as np

import Ofiura_Estimation as o_e
import Ofiura_planning as o_p
import Ofiura_ApriorPlanning as o_ap
import Ofiura_Qualitat as o_q

"""
Экстрагируем два параметра диода: N, коэфф. неидеальности и Is
"""

def diode (x,b,c=None):
    """
    Простая функция, выводящая ток через диод с параметрами в векторе b и при входном напряжении x[0]
    """
    Vd=x[0] #напряжение на диоде
    Is=b[0]
    N=b[1]
    #FT - температурный потенциал
    # g=1.60217657E-19 #Кл, заряд электрона
    # K=1.380648813E-23 #постоянная Больцмана
    # T=273+27 #температура p-n перехода в Кельвинах
    #
    #
    # FT=K*T/g #0.0258
    FT=0.02586419 #подогнанное по pspice
    I=Is*(math.exp(Vd/(FT*N)) -1)
    return [I]

def jacdiode (x,b,c=None,y=None):
    """
    Возвращает якобиан - по Is и N
    """
    FT=0.02586419 #подогнанное по pspice

    return np.matrix([math.exp(x[0]/(FT*b[1])) - 1, -b[0]*x[0]*math.exp(x[0]/(FT*b[1]))/(FT*b[1]**2) ]  )

def testDiode():

    b=[1e-14, 1]
    x=[10]
    print (jacdiode(x,b))


    print(diode(x=x,b=b))

    print('difference in percents:')
    pspice=8.19317e153 #10 V
    #pspice=12.9482 #0.9 V
    #pspice=9.05161E+069 #5 V

    print(100*(pspice-diode(x=x,b=b)[0])/pspice)

#0.0026638081177255196
    rng=np.arange(0.01,1,0.01)
    #снимем ВАХ
    resrng=[diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    plt.plot(rng , resrng)
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()

def testDiodeParameterExtraction():
    """
    пробуем экстрагировать коэффициенты из модели диода
    коэффициенты модели: Ток утечки Is, коэффициент неидеальности N, омическое сопротивление, параллельное диоду R
    входные параметры: напряжение, приложенное источником к системе резистор-диод
    +-----------|||||---------->|--------- -
    Резистор подключен до диода
    :return:
    """

    jacf=jacdiode
    funcf=diode

    #теперь попробуем сделать эксперимент.
    c={}
    Ve=np.asmatrix( [0.1]   )


    btrue=[1.238e-14, 1.8]
    bstart=np.array(btrue)-np.array(btrue)*0.2
    bend=np.array(btrue)+np.array(btrue)*0.2
    binit=[1e-10,1.1]


    xstart=[0.01]
    #xend=[20,60]
    xend=[2]

    N=30
    print("performing normal research:")
    startplan =  o_p.makeUniformExpPlan(xstart, xend, N)


    o_p.writePlanToFile(startplan)

    measdata = o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6, sign=1)
    #как мы помним, в случае неявных функций должно ставить sign=0
    print (gknu[0])
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))


    N=20
    print("performing aprior plan:")
    oplan=o_ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, funcf, Ntries=6)[1]
    o_p.writePlanToFile(oplan)
    measdata = o_p.makeMeasAccToPlan(funcf, oplan, btrue, c,Ve )
    gknu=o_e.grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6, sign=1)
    print (gknu[0])
    print (o_q.getQualitat(measdata, gknu[0], Ve,  funcf, c))




testDiodeParameterExtraction()