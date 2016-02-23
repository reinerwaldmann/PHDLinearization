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
from Cases.CasesUtilStuff import IterationInfoAcceptor
import Fianora.Fianora_Estimators as f_e


class DiodeMainScriptMantissaEx2(Fianora_MainScripts.DiodeMainScript):
    """
    Главный скрипт инкапсулирован в данном классе
        в init задаётся вычислительная задача, proceed может быть стандартным, может
        перопределяться

    """
    def __init__(self):

        Fianora_MainScripts.DiodeMainScript.__init__(self)

        estimator = f_e.NGEEstimatorUnlimited()
        #estimator = f_e.NGEstimator()
        estimator.init_parameters(self.ec, self.model)

        self.estimator  = f_e.ConsoleEstimatorDecorator(estimator)

    # Контекст (estimator context) пишется на этапе инита, опции оценки берутся в момент просида


import math
def norm1 (b1,b2):
    absdif = [math.fabs(b11-b22) for b11,b22 in zip(b1,b2)]
    return max(absdif)

def norm2 (b1,b2):
    return np.linalg.norm(np.array(b1)-np.array(b2))

def norm3 (b1,b2):
    absdif = [math.fabs(b11-b22) for b11,b22 in zip(b1,b2)]
    return sum(absdif)


def main():


    rs=[]
    nirc=205 # рассматриваемое число кейсов


    iiadumplist=[]


    for i in range (nirc):
        try:
            dms = DiodeMainScriptMantissaEx2()

            iia = IterationInfoAcceptor()

            dms.options.lst_data_abt_iteration = iia
            dms.proceed()
            # получили список значений b на каждой итерации

            btrue = dms.ec.btrue


            #Здесь определяется норма <<<>>>
            # такщем-та вроде как выводятся гистограммы для каждой итерации начиная с 4-й по эвклидовой норме
            # iia - это список данных по каждой итерации

            normfunc = norm1
            #list_of_errors = [normfunc(btrue, m['b']) for m in iia]

            # for k in iia:
            #     print (k, '\n')



            #plt.plot([z['b'][0] for z in iia])
            #plt.plot([z['b'][1] for z in iia])

            #plt.plot([z['b'][2] for z in iia])
            #print (iia)

            iiadumplist.append(iia)

            # fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
            # ax.plot([z['b'][1] for z in iia])
            # #это у нас графики значений компонента вектора по итерациям.
            # fig.savefig ('resfiles/M2/graphs/img_{0}.png'.format(i))
            # plt.close(fig)

            #получили список ошибок

            # print ('len', len(list_of_errors))
            # #if len(list_of_errors)<20: #работать только с теми случаями, где кол-во итераций не превысило обычное
            # #    #rs.append(list_of_errors)
            # print (list_of_errors)
            # rs.append(list_of_errors)

            #всадили его в результат
        except:

            pass


    # # складываем все результаты в пикулёвую БД в файл. Так-то б в мускул сложить, но то надо табличечку делать

    # # проще всего вытянуть из пикуля


    import pickle
    with open('resfiles/resdump205_DISP__.dat', 'wb') as f:
        pickle.dump(iiadumplist, f)





    # теперь rs это матрица список ошибок, прикол в том, шо каждая всунутая строка с инфой по итерациям
    # какбы чуть разной длины
    # берём срез в другую сторону
    rsT=[]
    for i in range (20): #свыше 20 итераций уже не рассматриваем, нах, до 3й тоже неинтересно
        try:
            rsT.append([rs[k][i] for k in range(nirc) ])
        except:
            break
    #теперь есть список списков погрешностей на первой, второй и далее итерациям
    #фигарим гистограммы за каждый список



    #[print (d) for d in rs]
    #print ()
    #[print (d) for d in rsT]

    ni=0


    import os

    filelist = [ f for f in os.listdir("resfiles/M2/hists") if f.endswith(".png") ]
    for f in filelist:
        try:
            os.remove("resfiles/M2/hists/"+f)
        except:
            pass

    for d in rsT:
        #d = [x for x in d if x<.1]
        ni+=1
        fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis

        d = list(filter( lambda x: x<.01, d))
        print (d)
        ax.hist(d,50)
        fig.savefig ('resfiles/M2/hists/img_{0}.png'.format(ni))
        plt.close(fig)






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

    # list_of_errors = [np.linalg.norm(np.array(btrue)-np.array(m['b'])) for m in iia]
    # print (list_of_errors)
    # plt.plot( list_of_errors)
    # plt.show()






 #   plt.show()



    # построить три графика, как падает погрешность





if __name__ == '__main__':
    main()



