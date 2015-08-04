"""
Экспериментальный вариант: ещё одна попытка автоматически делать производные
"""
__author__ = 'vasilev_is'
__date__ = "3 Aug 2015"
__version__ = "$Revision: 0 $"
__credits__ = "Nastya, Rammstein"


import math
import numpy as np
import sympy
import Fianora.Fianora_Models


class StringEMTransistorModel (EMTransistorModel):
    def evalfunc(self,dct, evalfn):
        """
        Рекурсивное вычисление словаря любой вложенности
        :param dct: строка или словарь, может быть словарь словарей и так далее
        :param evalfn: функция, которая принимает на вход строку, возвращает вычисленное значение
        :return: тот же словарь, но вместо строк нечто вычисленное
        """
        # при этом eval_string должна принимать на вход лишь строку

        if type(dct) is str:
            return (evalfn(dct))
        for k,v in dct.items():
            dct[k] = self.evalfunc(v, evalfn)

        return dct


    def make_deriv_components (self, complist, parameter_str_list):
        """
        Делает производные по компонентам - принимает на вход словарь компонентов и список параметров, возвращает словарь
        их производных
        По сути принимает на вход словарь хоть чего и у значений этого чего берёт производные

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


    def evalfn(self, _str, x, b, y):
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

        return eval(_str,  locals(), globals())     # исполняет полученную строку


    def model_from_components(self, **cval):
        """
        *
        Определяет, как из компонентов строится модель
        Даём на вход именованными аргументами значения вычисленных компонентов - на выходе получаем значение модели

        Принимает на вход компоненты в виде словаря, по факту как именованные аргументы, а выдаёт выходной вектор модели

        :param **cval: словарь значений компонентов
        :return: вектор значений модели. А может там и производной, мы-то не знаем =))
        """
        Ie = cval['Icc'] + cval['Ibe']
        Ic = cval['Ibc'] - cval['Icc']
        Ib = Ie-Ic

        return Ie, Ic, Ib


    def main_script_jac(x, b, y):
        """
        *
        Типа главный скрипт, который выбрасывает производную в зависимости от x, b, y
        В нём определены такие входные данные:
        список параметров
        определение компонентов через те параметры и не только
        """

        parameter_str_list = ['IS', 'BF', 'NR', 'NF', 'BR']     # список параметров



        #   компоненты
        Icc = 'IS * (exp(Vbe/(NF*FT))-exp(Vbc/(NR*FT) ))'
        Ibe = '(IS/BF) * (exp(Vbe/(NF*FT))-1)'
        Ibc = '(IS/BR) * (exp(Vbc/(NR*FT))-1)'

        # компоненты перевод в строку
        complist = {'Icc': Icc, 'Ibe': Ibe, 'Ibc':Ibc}


        #   общие параметры, которые можно определить по модели
        k = 3   # число уравнений
        m = 5   # число коэффициентов


        deriv_comp = make_deriv_components (complist, parameter_str_list)   # получили словарь производных по компонентам
        evfn = lambda str : evalfn(str, x, b, y)
        deriv_comp_estimated =  evalfunc(deriv_comp, evfn)  # взяли тот словарь и всё повычисляли

        res = np.matrix(k, m)   # строки столбцы

        i = 0
        for par in parameter_str_list:
            col = model_from_components(Icc = deriv_comp_estimated[par]['Icc'], Ibe = deriv_comp_estimated[par]['Ibe'], Ibc = deriv_comp_estimated[par]['Ibc'])
            for j in range(len(col)):
                res[j][i] = col[j]
        return res
















    #print (help(__name__))