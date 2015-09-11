"""
 2. Получение данных о точности результата на каждой итерации
        проведение 100 испытаний SimpleDiode составление таблицы
        номер итерации/число верных цифр для каждого b с сортировкой по возр, типа словарь/погрешность всего вектора/
        /длина мантиссы по погрешности с учётом порядка числа/Sk
"""

__author__ = 'vasilev_is'

import numpy as np
import matplotlib.pyplot as plt

from Fianora import Fianora_MainScripts



import Fianora.Fianora_Models as f_m
import Fianora.Fianora_Estimators as f_e
import Fianora.Fianora_Measurers as f_me
import Fianora.Fianora_Planners as f_p
import Fianora.Fianora_static_functions as f_sf


class IterationInfoAcceptor ():
    """
    Это класс, который ведёт себя по сути как список
    он принимает на вход любую хрень, и складирует её в список,
    тот список можно получить через переопределённые стандартные операторы

    """
    def __init__(self):
        self.lst_of_data=[]

    def accept (self, *args, **kwargs):
        if args and kwargs:
            self.lst_of_data.append((args,kwargs))
            return

        if args:
            self.lst_of_data.append(args)

        if kwargs:
            self.lst_of_data.append(kwargs)

    def __getitem__(self, key):
        return self.lst_of_data[key]

    def __iter__(self):
        return self.lst_of_data.__iter__()




class DiodeMainScriptMantissaEx2(Fianora_MainScripts.DiodeMainScript):
    """
    Главный скрипт инкапсулирован в данном классе
        в init задаётся вычислительная задача, proceed может быть стандартным, может
        перопределяться

    """
    def __init__(self):
        Fianora_MainScripts.DiodeMainScript.__init__()

        self.iia = IterationInfoAcceptor()
        self.options.lst_data_abt_iteration = self.iia



    def proceed(self, nocacheplan = False):
            return Fianora_MainScripts.DiodeMainScript.proceed(self,nocacheplan), self.iia

    # Контекст (estimator context) пишется на этапе инита, опции оценки берутся в момент просида






def main():
    pass




if __name__ == '__main__':
    main()



