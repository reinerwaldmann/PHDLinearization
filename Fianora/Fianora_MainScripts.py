
__author__ = 'vasilev_is'

import Fianora.Fianora_Models as f_m
import Fianora.Fianora_Estimators as f_e
import Fianora.Fianora_Measurers as f_me
import Fianora.Fianora_Planners as f_p
import Fianora.Fianora_static_functions as f_sf

import numpy as np



class AbstractMainScript():
    def __init__(self): # потому, что в конструкторе и делается самое интересное
        """
        Должен быть переопределён, ибо в нём конфигурируется вся радость
        :return:
        """
        pass

    def proceed(self, nocacheplan = False):
        # Вот здесь бы и фабричный метод пользовать, потому что получается, что этот методв
        # требует некоторых атрибутов. Явно видно, какие методы надо переопределить,
        # а вот насчёт атрибутов не всё так прозрачно.

            plan = self.planner.give_us_a_plan(nocacheplan)


            self.measdata = self.plan_measurer.measure_acc_to_plan(plan)
            self.estimator.estimate(self.measdata, f_e.Options())




class DiodeMainScript(AbstractMainScript):
    def __init__(self):

        btrue = [7.69e-8, 1.45 ,.0422] #номинальные значения диода D1N4001 с сайта, вроде официальной модели производителя
        binit = btrue
        Ve=np.array([[1.9e-5]])
        bstart=np.array(btrue)-np.array(btrue)*0.2
        bend=np.array(btrue)+np.array(btrue)*0.2
        xstart=[0.001]
        xend=[1]
        N=20
        ec =f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N)


        self.model = f_m.SimpleDiodeModel ('Diode_1N') # сделали модель
        self.measurer = f_me.ModelMeasurer(ec.Ve, self.model, ec.btrue) # сделали измерителя
        self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану
        #self.planner = f_p.DOptimalPlanner(self.N, self.xstart, self.xend, self.bstart, self.bend, self.model, verbose=1)
        self.planner = f_p.UniformPlanner(ec)

        self.measdata=None

        #формирование цепочки
        estimator = f_e.NGEstimator()
        estimator.init_parameters(ec, self.model)

        self.estimator  = f_e.ConsoleEstimatorDecorator(estimator)

class TransistorMainScript(AbstractMainScript):
     def __init__(self):
         #     .MODEL KT315 NPN (IS=10F BF=584.517 VAF=100 IKF=29.2714M ISE=131.803P
         # 3:  + NE=2.08337 BR=1.95214 IKR=9.99996M ISC=100.316P RE=0 RC=0 RB=0
         # 4:  + CJE=27.3893P VJE=700.001M MJE=500.287M CJC=27.3893P VJC=700.001M
         # 5:  + MJC=500.287M TF=450.287P XTF=499.984M VTF=10 ITF=10.2268M TR=153.383P)
         #
        inf=10e10
        IS = 10e-15
        BF = 584.517
        VAF = 100
        VAR = inf
        IKF =  29.2714e-3     # ток перехода к высококу уровню инжекции inf
        ISE = 131.803e-12      # генерационно-рекомбинационный ток насыщения  эмиттерного перех 0
        NE = 2.08337       # коэфф. неидеальности ген-рек тока эмиттерного перех   1
        NR = 1       # коэфф неидеальности для диффузного тока в инв режиме  1
        NF = 1       # коэфф неидеальности для диффузионного тока в нормальном режиме        1
        NC = 1       # коэфф неидеальности генерационно-рекомбинацоинного тока коллектора    1
        BR = 1.95214      # инверсный коэфф усиления тока в схеме с ОЭ 1
        IKR = 9.99996e-3     # ток перехода к высокому уровню инжекции в инверсном включении inf
        ISC = 100.316e-12     # генерационно-рекомбинационный ток насыщения колекторного перех 0
        RE = 1      # сопротивления эмиттера, коллектора, базы 0 0 0
        RC = 5.48635
        RB = 0


        parameter_str_list = ['IS', 'BF', 'NR', 'NF', 'BR']     # список параметров
        btrue = [IS, BF, NR, NF, BR]
        binit = btrue

        Ve = np.diag([1.9e-5]*3)



        bstart = np.array(btrue)-np.array(btrue)*0.2
        bend = np.array(btrue)+np.array(btrue)*0.2

        xstart = np.array([0.001, 0.001])
        xend = np.array([1, 1])
        N = 15

        ec = f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N) #упаковывание в контекст

        self.model = f_m.StringEMTransistorModel ('Q_NPN_KT513')

        ec.model = self.model




        self.measurer = f_me.ModelMeasurer(ec.Ve, self.model, ec.btrue) # сделали измерителя
        self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану


        self.planner = f_p.DOptimalPlanner(ec, 'cache', 1)


        #self.planner = f_p.UniformPlanner(ec)


        #формирование цепочки
        estimator = f_e.NGEstimator()
        estimator.init_parameters(ec, self.model)

        self.estimator  = f_e.ConsoleEstimatorDecorator(estimator)



def test():
    dm = TransistorMainScript()
    dm.proceed(nocacheplan=True)


test()