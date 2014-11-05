__author__ = 'reiner'
#попробуем промоделировать что попроще - диод. У него тоже нелинейная характеристика, и притом экспоненциальная
import math
import sympy
import Ofiura_planning as o_p
import numpy as np
import matplotlib.pyplot as plt

import scipy.optimize



def func (x,b,c=None):
    Is=b[0] #ток <strike>упячки</strike> утечки
    FT=0.026 #температурный потенциал
    Vin=x[0] #напряжение на входе

    return [Is*(math.exp(Vin/FT) - 1)]
def jac (x,b,c=None,y=None):
    Vin=x[0]
    FT=0.026

    return math.exp(Vin/FT)-1

def testNormal():
    #строим для проверки график
    rng=np.arange(0.01,2,0.01)
    resrng=[func([x],[1.28e-15])[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на колллекторе - снимаем ток коллектора.

    plt.plot(rng , resrng)
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()

#теперь последовательно добавляем резистор
#получается примерно вот что

def rSubfunc(Vel, V, R):
    #оптимизируемая функция - так мы найдём Vel
    return (V-Vel)/R-func([V],[1.28e-15])[0]

def funcWithR(x,b,c=None):
    #генеральная с оптимизатором
    V=x[0]  #Volts
    R=0.005 #Ohms

    ftoopt = lambda Vel: rSubfunc(Vel, V, R) #в данном случае оптимизируемая функция

    sol = scipy.optimize.minimize (ftoopt, [1])


    Vel = sol.x #ежели получилось, разумеется

    print("Vel=", Vel)
    print("Vr=", V-Vel)

    print ('Curr=',(V-Vel)/R)
    print ('Curr1=',func ([Vel], [1.28e-15])[0])

    return (V-Vel)/R #и возвращаем ток, наконец






rng=np.arange(2,3,0.01)
resrng=[funcWithR([x],[1.28e-15])[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на колллекторе - снимаем ток коллектора.

resNoR=[func([x],[1.28e-15])[0] for x in rng]

plt.plot(rng , resrng)
plt.plot(rng , resNoR)

#plt.axis([0.0,1.0,0,5])
plt.grid()
plt.show()

