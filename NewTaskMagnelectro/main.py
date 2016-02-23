"""
Файл содержит исследование оценки параметров некоего магнитоэлектрического коэффициента
модель наследуется напрямую от абстрактной. Производная определяется численно

Оценка параметров двумя вариантами: нашим, через Гаусса-Ньютона с ограничениями
и симплексным Нельдера-Мида, ибо не требует производной.


=Описание переменных=
==Параметры:

У Попова так:
B2(1) = 8.2 ' S11=8.2E-12 - эффективная податливость композита
B2(2) = -2.9 ' S12=-2.9E-12 -податливость композита
B2(3) = 6 'RO=6E3 - плотность
B2(4) = 5 'R=5E-3 - радиус диска
B2(5) = 9.2 'D31=9.2E-11- эффективный пьезоэлектрический (пьезомагнитный) модуль
B2(6) = 1.54 'E33=1.54E-8- диэлектрическая проницаемость
B2(7) = 6 'Q31=6E-11 - - пьезомагнитный модуль

v=-B2/B1

В нашем случае мы зададим B2 = V и получим

btrue=[     8.2e-12,
            0.3536585,
            6.0e3,
            5.0e-3,
            9.2e-11,
            1.54e-8,
            6.0e-11]

Ограничения:
b[i]>0 для всех

Коридор:
Для начала 30% в обе стороны с сохраняющимся ограничением на больше 0
Суммарный коридор 60%
(b[i]<btrue[i]+btrue[i]*.3)

==Входные переменные
w - что это такое, нет информации

Диапазон
346000 - 351000 шаг 0.5

==Дисперсия ошибок
V=[.0000001]
"""
from cmath import sqrt

import numpy as np


__author__ = 'reiner'



from Fianora.Fianora_Models import AbstractModel
import math


class Magnitoelectro_model (AbstractModel):
    def __init__(self):
        self.name = "Magnitoelectro. bessel: scipy approx. deriv: scipy approx (numerical)"

    def funcf(self, x,b):
        """
        Main model function
        returns y vector
        """
        from scipy.special import jv
        w=x[0]
        k2=sqrt(b[2]*b[0]*(1-b[1]**2)).real*w*2*math.pi*b[3]

        j0 = jv(0,k2)
        j1 = jv(1,k2)

        kp2=2*b[4]**2/(b[5]*b[0]*(1-b[1]))
        dr=k2*j0-(1-b[1])*j1
        da=1-kp2+kp2*(1+b[1])*j1/dr
        ael = (1/da) * (2*b[4]*b[6]/(b[5]*b[0]*(1-b[1]))) * ((dr-(1+b[1])*j1)/dr)

        return [1/ael]

    def jacf(self, x, b, y=None):
        #here comes a guide to write a numerical derivative function
        #
        # http://docs.scipy.org/doc/scipy-0.16.1/reference/generated/scipy.misc.derivative.html
        return None

    #
    #     def jacf(self, x, b, y=None):
    #         """
    #         Jacobian
    #         """
    #         raise NotImplementedError("jacf Please Implement this method")
    #


def getXYfile ():
    btrue=[8.2e-12,
            0.3536585,
            6.0e3,
            5.0e-3,
            9.2e-11,
            1.54e-8,
            6.0e-11]

    #x[0]=346000 - 351000 шаг 0.5
    x = list(np.arange(346000.0,351000.0,.5))
    y = list()
    memodel = Magnitoelectro_model  ()

    with open ("../resfiles/NewTaskMagnelectroXYVals.txt", "wt") as out_file:

        for xx in x:
            print ([xx],  memodel.funcf([xx], btrue), file=out_file)


            #y.append(memodel.funcf(x, btrue))



