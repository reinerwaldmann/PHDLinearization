from abc import abstractmethod
from abc import ABCMeta

__author__ = 'vasilev_is'
import math
from scipy import optimize
import numpy as np
import sympy

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
            solx=optimize.root(self.mf, solvinit, args=(x,b), jac=self.dfdy, method='lm').x
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


        IS = b[0]   # сквозной ток насыщения   1e-16
        BF = b[1]   # максимальное значение нормального коэфф усиления по току с ОЭ 100
        NR = b[2]       # коэфф неидеальности для диффузного тока в инв режиме  1
        NF = b[3]       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
        BR = b[4]      # инверсный коэфф усиления тока в схеме с ОЭ 1

        #
        # IS = b[0]       # сквозной ток насыщения   1e-16
        # BF = b[1]       # максимальное значение нормального коэфф усиления по току с ОЭ 100
        # VAF = b[2]      # прямое напряжение Эрли    inf
        # VAR = b[3]      # инверсное напряжение Эрли inf
        # IKF =  b[4]     # ток перехода к высококу уровню инжекции inf
        # ISE = b[5]      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
        # NE = b[6]       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
        # NR = b[7]       # коэфф неидеальности для диффузного тока в инв режиме  1
        # NF = b[8]       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
        # NC = b[9]       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
        # BR = b[10]      # инверсный коэфф усиления тока в схеме с ОЭ 1
        # IKR = b[11]     # ток перехода к высокому уровню инжекции в инверсном включении inf
        # ISC = b[12]     # генерационно-рекомбинационный ток насыщения колекторного перех 0
        # RE = b[13]      # сопротивление эмиттера
        # RC = b[14]      # сопротивление коллектора
        # RB = b[15]      # сопротивление базы


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

    def jacf(self, x, b, y=None):
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


class StringDefinedModel(AbstractModel):
    """
    Это абстрактный класс модели, которая задаётся строками
    Ключевая идея - модель задаётся как некоторая комбинация из компонентов.
    Вот методы, которые надо переопределять:

        make_parameter_str_list(self): - должен вернуть список строк - список параметров в том виде, в котором они
        будут встречаться в модели и в той последовательности, в которой они в векторе

        def make_modelstr(self):




    """


    def __init__(self, name):
        #паттерн "Фабричный метод" во все поля!

        AbstractModel.__init__(self, name)

        self.parameter_str_list = self.make_parameter_str_list()
        self.modelstr = self.make_modelstr()
        self.components = self.make_components()
        self.m = len(self.parameter_str_list)
        self.k = len(self.modelstr)  # число уравнений

    def make_parameter_str_list(self):
        raise NotImplementedError('pls implement me')

    def make_modelstr(self):
        raise NotImplementedError('pls implement me')

    def make_components(self):
        raise NotImplementedError('pls implement me')

    def eval_string(self, _str, x, b, y=None):
        """
        *
        Вычисляет значение строки, пользуясь тем, что знает, какой переменной соответствует
        тот или иной компонент вектора
        :param _str: строка
        :param x: вектор x
        :param b: вектор b
        :return:
        """
        raise NotImplementedError('pls implement me')


    def model_from_components (self, cval):
        """
        возвращает список вычисленных значений в условиях всех локальных переменных (т. е. определённых компонентов),
         и глобальных (в смысле подключеных модулей, главным образом)

        :param cval: словарь компонентов компонент - значение
        :return:
        """

        return [eval(st, cval, globals()) for st in self.modelstr]

    def evalfunc(self, dct, evalfn):
        """
        Рекурсивное вычисление словаря любой вложенности
        :param dct: строка или словарь, может быть словарь словарей и так далее
        :param evalfn: функция, которая принимает на вход строку, возвращает вычисленное значение
        :return: тот же словарь, но вместо строк нечто вычисленное
        """
        # при этом eval_string должна принимать на вход лишь строку


        if type(dct) is str:
            return (evalfn(dct))

        if type(dct) in [float, int]:
            return dct

        res={}
        for k,v in dct.items():
            res[k] = self.evalfunc(v, evalfn)

        return res

    def make_deriv_components (self, complist, parameter_str_list):
        """
        Делает производные по компонентам - принимает на вход словарь компонентов и список параметров, возвращает словарь
        их производных
        По сути принимает на вход словарь хоть чего и у значений этого чего берёт производные
        Не обязательно вообще делить модель на компоненты  -  можно и у уравнений сразу брать производные

        :param complist: словарь компонентов - идентификатор компонента - его строковое представление
        :param parameter_str_list: список параметров, строковый
        :return: словарь параметр - компонент - производная в строковом виде
        """

        rs = {}
        for par in parameter_str_list:
            # print()
            # print(par)
            cd = {}
            for name, comp in complist.items():
                # pre-processing of string
                comp = comp.replace("math.exp", "exp")

                sd = sympy.diff(comp, par)
                sd = str(sd)

                # post-processing of string
                sd = sd.replace("exp", "math.exp")

                cd[name] = sd
            rs[par] = cd

        return rs

    def jacf(self, x, b, y=None):
        """
        Типа главный скрипт, который выбрасывает производную в зависимости от x, b, y
        """

        deriv_comp = self.make_deriv_components (self.components, self.parameter_str_list)   # получили словарь производных по компонентам
        evfn = lambda str : self.eval_string(str, x, b, y)
        deriv_comp_estimated = self.evalfunc(deriv_comp, evfn)  # взяли тот словарь и всё повычисляли




        res = np.zeros((self.k, self.m))   # строки столбцы

        i = 0
        for par in self.parameter_str_list:
            col = self.model_from_components(deriv_comp_estimated[par])
            for j in range(len(col)):
                res[j][i] = col[j]
            i += 1
        return res

    def funcf(self, x,b):
        evfn = lambda _str: self.eval_string(_str, x, b)

        components_estimated = self.evalfunc(self.components, evfn)
        return np.array(self.model_from_components(components_estimated))


class StringEMTransistorModel (SemiconductorModel, StringDefinedModel):
    """
    Это попытка создать модель на основе строковой записи компонентов
     и определения, как из этих компонентов должна создаваться модель

     Определяемые подклассом сущности:
     +список строковых параметров
     +строковое представление того, как из компонентов строится модель в виде вектора уравнений, в которые входят компоненты
     +словарь компонентов. Компоненты уже могут определяться через числа, константы и так далее
     +функция, которая со значениями векторов x,b,y вычисляет значение поданной строки


    """
    #def __init__(self, name):
        #EMTransistorModel.__init__(name)
        #StringDefinedModel.__init__(name)
        #AbstractModel.__init__(name)



    def make_parameter_str_list(self):
        return ['IS', 'BF', 'NR', 'NF', 'BR']     # список параметров

    def make_modelstr(self):
        return  ['Icc+Ibe',
                         'Ibc-Icc',
                         'Icc+Ibe-Ibc+Icc']   # строковое представление того, как из компонентов строится модель

    def make_components(self):
        return {'Icc': 'IS * (math.exp(Vbe/(NF*FT))-math.exp(Vbc/(NR*FT) ))',
                            'Ibe': '(IS/BF) * (math.exp(Vbe/(NF*FT))-1)',
                            'Ibc': '(IS/BR) * (math.exp(Vbc/(NR*FT))-1)'}

    def eval_string(self, _str, x, b, y=None):
        """
        *
        Вычисляет значение строки, пользуясь тем, что знает, какой переменной соответствует
        тот или иной компонент вектора

        :param _str: строка
        :param x: вектор x
        :param b: вектор b
        :return:
        """


        b = list(map(float,b))
        x = list(map(float,x))
        Vbe = -1*x[0]
        Vbc = -1*x[1]
        IS = b[0]   # сквозной ток насыщения   1e-16
        BF = b[1]   # максимальное значение нормального коэфф усиления по току с ОЭ 100
        NR = b[2]       # коэфф неидеальности для диффузного тока в инв режиме  1
        NF = b[3]       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
        BR = b[4]      # инверсный коэфф усиления тока в схеме с ОЭ 1
        FT = self.FT

        ev = eval(_str,  locals(), globals())


        return ev     # исполняет полученную строку




