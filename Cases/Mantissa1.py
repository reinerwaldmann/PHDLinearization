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
    def __init__(self, filename=None):
        self.lst_of_data=[]
        if filename:
            self.filename = filename
            with open(self.filename, 'w') as f:
                pass


    def accept (self, *args, **kwargs):
        if self.filename:
                with open(self.filename, 'a') as f:
                    f.write(args.__str__()+"\t"+kwargs.__str__()+"\n")



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


    def proceed(self, nocacheplan = False):
        return Fianora_MainScripts.DiodeMainScript.proceed(self,nocacheplan), self.iia

    # Контекст (estimator context) пишется на этапе инита, опции оценки берутся в момент просида




def main():
    dms = Fianora_MainScripts.DiodeMainScript()
    iia = IterationInfoAcceptor('resfiles/Mantissa2res.txt')
    dms.options.lst_data_abt_iteration = iia
    dms.proceed()
    btrue = dms.ec.btrue
    print (btrue)


    #получает графики убывания погрешностей в зависимости от номера итерации по
    # теперь получить таблицу погрешностей абсолютных по итерациям
    list_of_errors = [np.fabs(np.array(btrue)-np.array(m['b'])) for m in iia]

    flatlel = [ [lll[i] for lll in list_of_errors] for i in range(len(btrue))]

    ar = list(range(len(flatlel[0])))
    print (list_of_errors)

    for i in range (len(btrue)):
        plt.plot(ar,flatlel[i])
        plt.savefig ('resfiles/M2/{0}.png'.format(i))
        plt.show()






 #   plt.show()



    # построить три графика, как падает погрешность





if __name__ == '__main__':
    main()



