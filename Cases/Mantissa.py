"""
Серия экспериментов про мантиссу

  =Группа экспериментов Mantissa тред=
    Назначение: обеспечить практическую поддержку разрабатываемого раздела по мантиссе
    Эксперименты:
        1. Проверка гипотезы постоянства числа итераций
        берётся стандартный кейс (диод), число проверок (100) на каждой проверке случайно выбирается binit
        и случайно выбирается btrue, пишем в список число итераций, проверяем гипотезу распределения,
        строим диаграмму рассеяния, определяем среднее, дисперсию
        2. Получение данных о точности результата на каждой итерации
        проведение 100 испытаний SimpleDiode составление таблицы
        номер итерации/число верных цифр для каждого b с сортировкой по возр, типа словарь/погрешность всего вектора/
        /длина мантиссы по погрешности с учётом порядка числа/Sk
        3.MPM DebugMeasurer MPM,
        DiodeModelAdvancedMPM
        NGEstimatorMPM с поддержкой переменной мантиссы
        проверять GMP как основу MPM
        Прогон с большой мантиссой расширенного диода, прогон с переменной мантиссой, выявление экономии времени
        и точностных характеристик


        Ахтунг!
        Дело в том, что при ограниченном вычислении есть внешнее число итераций, число итераций обёртки
        и внутреннее число итераций, то есть число итераций метода
        Оно бывает вообще на последней итерации какбе 0, но вот с диодом покамест получается что
        пару раз подряд крутит по 4-5 итераций.
        Так вот, перемена мантиссы и прочая дичь нужны именно для внутреннего числа итераций итерационного метода.
        Можно вообще блокировать внешний итерационный процесс

        Ещё ахтунг: нигде не обнаружено описание структуры данных выхлопа оценочной процедуры

"""

__author__ = 'vasilev_is'

import numpy as np
import matplotlib.pyplot as plt

from Fianora import Fianora_MainScripts


from Fianora import *

import Fianora.Fianora_Models as f_m
import Fianora.Fianora_Estimators as f_e
import Fianora.Fianora_Measurers as f_me
import Fianora.Fianora_Planners as f_p
import Fianora.Fianora_static_functions as f_sf


# 1

class DiodeMainScriptMantissaEx1(Fianora_MainScripts.AbstractMainScript):
    """
    Главный скрипт инкапсулирован в данном классе
        в init задаётся вычислительная задача, proceed может быть стандартным, может
        перопределяться

    """
    def __init__(self):
        Fianora_MainScripts.AbstractMainScript.__init__(self)
        self.options.verbose = 0
        self.options.verbose_wrapper = 0

        btrue = [7.69e-8, 1.45 ,.0422] #номинальные значения диода D1N4001 с сайта, вроде официальной модели производителя
        binit = btrue
        Ve=np.array([[1.9e-5]])
        bstart=np.array(btrue)-np.array(btrue)*0.4
        bend=np.array(btrue)+np.array(btrue)*0.4
        xstart=[0.001]
        xend=[1]
        N=100
        self.ec =f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N)

        self.model = f_m.SimpleDiodeModel ('Diode_1N') # сделали модель

        #self.planner = f_p.DOptimalPlanner(self.N, self.xstart, self.xend, self.bstart, self.bend, self.model, verbose=1)
        self.planner = f_p.UniformPlanner(self.ec)
        self.measdata=None
        #формирование цепочки
        self.estimator = f_e.NGEstimator()
        self.plan = self.planner.give_us_a_plan(nocache = True) # чтоб каждый раз не менять план, как оно делается,
        # когда это запускается в proceed



    def proceed(self, nocacheplan = False):
           # print (self.ec.__dict__)
            self.estimator.init_parameters(self.ec, self.model)
            self.measurer = f_me.ModelMeasurer(self.ec.Ve, self.model, self.ec.btrue) # сделали измерителя
            self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану
            self.measdata = self.plan_measurer.measure_acc_to_plan(self.plan)
            return self.estimator.estimate(self.measdata, self.options)






class TransistorMainScriptMantissaEx1(DiodeMainScriptMantissaEx1):
     def __init__(self):
         #     .MODEL KT315 NPN (IS=10F BF=584.517 VAF=100 IKF=29.2714M ISE=131.803P
         # 3:  + NE=2.08337 BR=1.95214 IKR=9.99996M ISC=100.316P RE=0 RC=0 RB=0
         # 4:  + CJE=27.3893P VJE=700.001M MJE=500.287M CJC=27.3893P VJC=700.001M
         # 5:  + MJC=500.287M TF=450.287P XTF=499.984M VTF=10 ITF=10.2268M TR=153.383P)
         #
        Fianora_MainScripts.AbstractMainScript.__init__(self)
        self.options.verbose = 1
        self.options.verbose_wrapper = 1

        inf=10e10
        #числовые значения изменены шоб было по 7 значимых
        IS = 11.23456e-15
        BF = 584.5171
        VAF = 112.3456
        VAR = inf
        IKF =  29.27141e-3     # ток перехода к высококу уровню инжекции inf
        ISE = 131.8031e-12      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
        NE = 2.083371       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
        NR = 1       # коэфф неидеальности для диффузного тока в инв режиме  1
        NF = 1       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
        NC = 1       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
        BR = 1.952141      # инверсный коэфф усиления тока в схеме с ОЭ 1
        IKR = 9.999961e-3     # ток перехода к высокому уровню инжекции в инверсном включении inf
        ISC = 100.3161e-12     # генерационно-рекомбинационный ток насыщения колекторного перех 0
        RE = 1      # сопротивления эмиттера, коллектора, базы 0 0 0
        RC = 5.48635
        RB = 0
        parameter_str_list = ['IS', 'BF', 'NR', 'NF', 'BR']     # список параметров
        btrue = [IS, BF, NR, NF, BR]
        Ve = np.diag([1.9e-5]*3)
        bstart = np.array(btrue)-np.array(btrue)*0.3
        bend = np.array(btrue)+np.array(btrue)*0.3
        binit = f_sf.uniformVector(bstart, bend)
        xstart = np.array([0.001, 0.001])
        xend = np.array([1, 1])
        N = 30

        self.ec = f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N) #упаковывание в контекст
        self.model = f_m.StringEMTransistorModel ('Q_NPN_KT513')



        self.measurer = f_me.ModelMeasurer(self.ec.Ve, self.model, self.ec.btrue) # сделали измерителя
        self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану
        self.planner = f_p.UniformPlanner(self.ec)

        #формирование цепочки
        self.estimator = f_e.NGEstimator()

        self.measdata=None

        self.plan = self.planner.give_us_a_plan(nocache = True) # чтоб каждый раз не менять план, как оно делается,




def main ():
    """
    Главная функция, осуществляющая эксперимент
    :return:
    """
     #bstart = [  6.15200000e-08  , 1.16000000e+00,   3.37600000e-02] # 20%
    #bend =   [  9.22800000e-08 ,  1.74000000e+00 ,  5.06400000e-02] # 20%


#    bstart = [  4.61400000e-08,   8.70000000e-01,   2.53200000e-02]  # 40%
 #   bend = [  1.07660000e-07,   2.03000000e+00,   5.90800000e-02]  # 40%

    # transistor

    inf=10e10
    #числовые значения изменены шоб было по 7 значимых

    IS = 11.23456e-15
    BF = 584.5171
    VAF = 112.3456
    VAR = inf
    IKF =  29.27141e-3     # ток перехода к высококу уровню инжекции inf
    ISE = 131.8031e-12      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
    NE = 2.083371       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
    NR = 1       # коэфф неидеальности для диффузного тока в инв режиме  1
    NF = 1       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
    NC = 1       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
    BR = 1.952141      # инверсный коэфф усиления тока в схеме с ОЭ 1
    IKR = 9.999961e-3     # ток перехода к высокому уровню инжекции в инверсном включении inf
    ISC = 100.3161e-12     # генерационно-рекомбинационный ток насыщения колекторного перех 0
    RE = 1      # сопротивления эмиттера, коллектора, базы 0 0 0
    RC = 5.48635
    RB = 0
    btrue = [IS, BF, NR, NF, BR]
    bstart = np.array(btrue)-np.array(btrue)*0.3
    bend = np.array(btrue)+np.array(btrue)*0.3

    from Cases.Mantissa1 import IterationInfoAcceptor

    iia = IterationInfoAcceptor ('resfiles/M1.txt')

    msa = TransistorMainScriptMantissaEx1()

    for i in range (50):

        btrue1 = f_sf.rangomNormalvariateVector(bstart, bend)

        msa.ec.btrue = btrue1
        msa.ec.bstart = bstart
        msa.ec.bend = bend

        try:
            print ('entering')
            rs = msa.proceed()
            iia.accept(numiter=rs['numiter'], Sk=rs['Sk'], AvDif=rs['AvDif'])


        except:
            print ('uno problemo')
            pass

        plt.hist ([x[0] for x in iia],20)
        plt.show()

        plt.savefig('resfiles/MantissaResTransistor.png')


if __name__=='__main__':
      main()




# l = [14, 14, 10, 27, 10, 8, 10, 11, 51, 10, 48, 19, 8, 5, 31, 15, 9, 10, 13, 5, 15, 15, 51, 34, 49, 6, 5, 10, 15, 46, 11, 10, 16, 22, 7, 11, 42, 7, 10, 29, 27, 15, 9, 7, 38, 27, 10, 29, 15, 12, 17, 7, 11, 8, 32, 51, 12, 23, 40, 10, 24, 34, 34, 40, 9, 7, 26, 51, 51, 51, 23, 14, 29, 22, 10, 50, 33, 8, 6, 6, 19, 11, 22, 9, 14, 30, 51, 14, 15, 13, 10, 6, 7, 17, 17, 23, 51, 17, 7, 10, 5, 7, 24, 29, 7, 21, 47, 7, 16, 20, 19, 23, 5, 22, 9, 7, 12, 9, 15, 32, 29, 13, 11, 16, 9, 16, 51, 20, 5, 8, 24, 46, 9, 25, 6, 10, 7, 36, 21, 10, 6, 11, 17, 15, 51, 8, 10, 11, 19, 30, 9, 51, 34, 9, 27, 22, 15, 14, 9, 45, 16, 29, 9, 10, 17, 10, 20, 51, 20, 17, 51, 11, 30, 20, 8, 11, 19, 41, 14, 12, 49, 14, 10, 5, 13, 8, 10, 28, 7, 7, 6, 8, 11, 8, 10]
#binit const btrue unif

#ll=[18, 7, 13, 13, 12, 9, 7, 51, 6, 9, 9, 5, 7, 6, 7, 8, 14, 7, 12, 13, 6, 12, 9, 14, 20, 18, 18, 7, 43, 22, 11, 29, 6, 10, 34, 9, 9, 13, 7, 5, 12, 10, 7, 15, 13, 9, 12, 6, 19, 6, 12, 33, 13, 9, 13, 5, 11, 9, 6, 7, 12, 24, 8, 11, 12, 7, 51, 5, 10, 6, 10, 17, 12, 18, 12, 32, 9, 13, 7, 8, 15, 7, 11, 30, 8, 23, 16, 7, 9, 6, 6, 10, 6, 8, 16, 6, 12, 16, 5, 9, 6, 8, 26, 8, 8, 6, 5, 8, 10, 9, 12, 33, 7, 8, 10, 9, 9, 17, 7, 25, 10, 35, 7, 13, 12, 7, 9, 6, 14, 4, 10, 16, 8, 8, 22, 16, 9, 7, 21, 22, 12, 11, 12, 7, 14, 8, 10, 17, 7, 7, 6, 16, 17, 23, 9, 29, 15, 6, 20, 7, 11, 37, 47, 20, 8, 8, 8, 13, 12, 29, 9, 25, 5, 8, 9, 11, 8, 17, 7, 20, 10, 5, 10, 12, 10, 18, 14, 22, 18, 12, 15, 9, 34, 11, 7, 8, 7, 20, 24, 6]
#binit const btrue norm

#plt.hist (ll,50)
#plt.show()



 # def hist(x, bins=10, range=None, normed=False, weights=None, cumulative=False,
#          bottom=None, histtype='bar', align='mid', orientation='vertical',
#          rwidth=None, log=False, color=None, label=None, stacked=False,
#          hold=None, **kwargs):

# import random
# rr = [random.normalvariate (1,100) for i in range(1000)]
# plt.hist (rr, 20)
# plt.show()


