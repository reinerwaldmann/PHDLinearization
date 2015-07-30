__author__ = 'vasilev_is'
import math
from scipy import optimize
import numpy as np

# абстракции из данного раздела можно раскрутить и для автоматизированного нахождения производных символьным методом



class AbstractModel:
    """
    Defines a model abstraction.
    """

    def __init__(self, name):
        self.name = name

    def funcf (self, x,b):
        """
        Main model function
        """
        raise NotImplementedError("funcf Please Implement this method")

    def jacf(self, x, b, y=None):
        """
        Jacobian
        """
        raise NotImplementedError("jacf Please Implement this method")

    def cal_funcf(self):
        return self.funcf

    def cal_jacf(self):
        return self.jacf




class SemiconductorModel(AbstractModel):
    FT = 0.02586419


class ImplicitModel(AbstractModel):
    """
    Defines model class, which model defined implicitly, defining 2 abstract functions: mf, which
    is model and dfdy and dfdb.
    Then it implements solver, and jacobian, as template methods (implemented method which use abstract methods)
    """

    def dfdy(self, y, x, b):
        raise NotImplementedError("hess Please Implement this method")

    def mf(self, y, x, b):
        raise NotImplementedError("hess Please Implement this method")

    def dfdb(self, y, x, b):
        raise NotImplementedError("hess Please Implement this method")

    def funcf(self, x, b):

        solvinit=[1]

        try:
            solx=optimize.root(self.mf, solvinit, args=(x,b,c), jac=self.dfdy, method='lm').x
            #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
            #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
        except BaseException as e:
            print ('solving function: Error in findroot=',e)
            return None

        if solx-solvinit==[0]*len(solx):
            return None

        return solx

    def jacf(self, x, b, y=None):

        # если значение y не приходит, то мы его находим вызовом солвера
        # кстати эта часть одинакова, с другой стороны,
        # в тех функциях, где модель задаётся явно, y не участвует в производной и потому можно всегда передавать None
        if not y:
            y = self.funcf(x, b)
        return np.dot(np.linalg.inv(self.dfdy(y, x, b)), self.dfdb(y, x, b))



class SimpleDiodeModel(SemiconductorModel, ImplicitModel): # множественное наследование - от полупроводниковой модели и от неявной модели
    def mf (self, y,x,b):
        # модельная функция
        mm=float(b[0]*(math.exp((x[0]-y[0]*b[2])/(self.FT*b[1])) -1)-y[0])
        return [mm]

    def dfdy(self,y,x,b):
        return np.array ([[ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(self.FT*b[1]))/(self.FT*b[1])]])

    def dfdb(self,y,x,b):
        return np.matrix( [ [math.exp((-b[2]*y[0] + x[0])/(self.FT*b[1])) - 1,
                                          -b[0]*(-b[2]*y[0] + x[0])*math.exp((-b[2]*y[0] + x[0])/(self.FT*b[1]))/(self.FT*b[1]**2),
                                          -b[0]*y[0]*math.exp((-b[2]*y[0] + x[0])/(self.FT*b[1]))/(self.FT*b[1])]])



class AdvancedDiodeModel(SemiconductorModel, ImplicitModel):
    pass


class EMTransistorModel(SemiconductorModel):
    def funcf(self, x,b):
        """
        Модель Эберса-Молла явная
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
        RE = b[13]      # сопротивление эмиттера
        RC = b[14]      # сопротивление коллектора
        RB = b[15]      # сопротивление базы


        FT = self.FT

        Icc = IS * (math.exp(Vbe/(NF*FT))-math.exp(Vbc/(NR*FT) ))
        Ibe = (IS/BF) * (math.exp(Vbe/(NF*FT))-1)
        Ibc = (IS/BR) * (math.exp(Vbc/(NR*FT))-1)
        # генерационно-рекомбинационные составляющие

        Ie = Icc+Ibe
        Ic = Ibc-Icc
        Ib = Ie-Ic

        #Ie = IS * (math.exp(Vbe/(NF*FT))-math.exp(Vbc/(NR*FT) )) +  (IS/BF) * (math.exp(Vbe/(NF*FT))-1)
        #Ic = IS * (math.exp(Vbe/(NF*FT))-math.exp(Vbc/(NR*FT) )) + (IS/BR) * (math.exp(Vbc/(NR*FT))-1)
        #ток базы нам не требуется, хотя..

        y=[Ie, Ic, Ib]
        return y




class GPTransistorModel(SemiconductorModel, ImplicitModel):
    def funcf(self, x, b):
        """
        Модель Гуммеля-Пуна
        :param x:
        :param b:
        :return:
        """
        import mpmath as mpm

        FT=self.FT

        #pr = lambda  x: mpm.mpf(str(x))
        pr=float

        b = list(map(pr  ,b))
        x = list(map(pr ,x))


        Vbe = -1*x[0]
        Vbc = -1*x[1]
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
