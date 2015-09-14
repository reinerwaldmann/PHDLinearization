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
    Если задан файл, то он будет туда писать по мере прихода данных
    При начальном обявлении файл перезаписывается

    """
    def __init__(self, filename=None, verbose=False):
        self.lst_of_data=[]
        self.verbose = verbose
        self.lineCounter = 0
        if filename:
            self.filename = filename
            with open(self.filename, 'w') as f:
                pass


    def accept (self, *args, **kwargs):
        if self.filename:
                with open(self.filename, 'a') as f:
                    f.write(self.lineCounter+"\t"+args.__str__()+"\t"+kwargs.__str__()+"\n")


        if self.verbose:
            print (self.lineCounter, args, kwargs)

        self.lineCounter+=1


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

        Fianora_MainScripts.DiodeMainScript.__init__(self)

        estimator = f_e.NGEEstimatorUnlimited()
        estimator.init_parameters(self.ec, self.model)

        self.estimator  = f_e.ConsoleEstimatorDecorator(estimator)

    # Контекст (estimator context) пишется на этапе инита, опции оценки берутся в момент просида




def main():
    #dms = Fianora_MainScripts.DiodeMainScript()
    dms = DiodeMainScriptMantissaEx2()
    iia = IterationInfoAcceptor('resfiles/Mantissa2res.txt')
    dms.options.lst_data_abt_iteration = iia
    dms.proceed()
    btrue = dms.ec.btrue
    print (btrue)




    #получает графики убывания погрешностей в зависимости от номера итерации ПО КОМПОНЕНТАМ
    # теперь получить таблицу погрешностей абсолютных по итерациям
    # list_of_errors = [np.fabs(np.array(btrue)-np.array(m['b'])) for m in iia]
    #
    # flatlel = [ [lll[i] for lll in list_of_errors] for i in range(len(btrue))]
    #

    # print (list_of_errors)
    #
    # for i in range (len(btrue)):
    #     plt.plot(flatlel[i])
    #     plt.savefig ('resfiles/M2/{0}.png'.format(i))
    #     plt.show()
    #

    #График убывания погрешности в зависимости В ЦЕЛОМ ПО НОРМАЛИ ВЕКТОРА

    list_of_errors = [np.linalg.norm(np.array(btrue)-np.array(m['b'])) for m in iia]
    print (list_of_errors)
    plt.plot( list_of_errors)
    plt.show()






 #   plt.show()



    # построить три графика, как падает погрешность





if __name__ == '__main__':
    main()



