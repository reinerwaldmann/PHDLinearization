
"""
Describes abstract main script
"""
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


        self.options = f_e.Options()
        self.options.verbose = 1
        self.options.verbose_wrapper = 1

        self.options.NSIG = 5
        self.options.NSIGGENERAL = 5

        self.ec = None




    def proceed(self, nocacheplan = False):
        # Вот здесь бы и фабричный метод пользовать, потому что получается, что этот методв
        # требует некоторых атрибутов. Явно видно, какие методы надо переопределить,
        # а вот насчёт атрибутов не всё так прозрачно.
            plan = self.planner.give_us_a_plan(nocacheplan)
            self.measdata = self.plan_measurer.measure_acc_to_plan(plan)
            return self.estimator.estimate(self.measdata, self.options)


class DiodeMainScript(AbstractMainScript):
    def __init__(self):

        AbstractMainScript.__init__(self)

        _btrue = [7.69e-8, 1.45 ,.0422] #номинальные значения диода D1N4001 с сайта, вроде официальной модели производителя
        Ve=np.array([[1.9e-20]])
        bstart=np.array(_btrue)-np.array(_btrue)*0.2
        bend=np.array(_btrue)+np.array(_btrue)*0.2

        binit=_btrue

        btrue = f_sf.rangomNormalvariateVector(bstart, bend)
        print (btrue)
        #btrue = _btrue




        xstart=[0.001]
        xend=[1]
        N=100
        self.ec =f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N)


        self.model = f_m.SimpleDiodeModel ('Diode_1N') # сделали модель
        self.ec.model = self.model


        self.measurer = f_me.ModelMeasurer(self.ec.Ve, self.model, self.ec.btrue) # сделали измерителя
        self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану

        #self.planner = f_p.DOptimalPlanner(self.ec)
        self.planner = f_p.UniformPlanner(self.ec)


        #self, ec, plancachefoldername='cache', verbose = False)



        self.measdata=None

        #формирование цепочки
        estimator = f_e.NGEstimator()
        estimator.init_parameters(self.ec, self.model)

        self.estimator  = f_e.ConsoleEstimatorDecorator(estimator)

class TransistorMainScript(AbstractMainScript):
     def __init__(self):
        AbstractMainScript.__init__(self)
         #     .MODEL KT315 NPN (IS=10F BF=584.517 VAF=100 IKF=29.2714M ISE=131.803P
         # 3:  + NE=2.08337 BR=1.95214 IKR=9.99996M ISC=100.316P RE=0 RC=0 RB=0
         # 4:  + CJE=27.3893P VJE=700.001M MJE=500.287M CJC=27.3893P VJC=700.001M
         # 5:  + MJC=500.287M TF=450.287P XTF=499.984M VTF=10 ITF=10.2268M TR=153.383P)
         #
        inf=10e10


        #IS = 10e-15
        IS = 1
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

        b_nominal = [IS, BF, NR, NF, BR]
        bstart = np.array(b_nominal)-np.array(b_nominal)*0.3
        bend = np.array(b_nominal)+np.array(b_nominal)*0.3

        #btrue = f_sf.rangomNormalvariateVector(bstart, bend)
        btrue =  [1.1204058410408533, 587.71110706589764, 1.0110737035569517, 0.97217262436055318, 1.7709083403162802]


        #btrue  = b_nominal
        binit = b_nominal

        print ('Btrue set to:', btrue)


        Ve = np.diag([1.9e-5]*3)


        xstart = np.array([0.001, 0.001])
        xend = np.array([.6, .6])
        N = 60

        ec = f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N) #упаковывание в контекст

        self.model = f_m.StringEMTransistorModel ('Q_NPN_KT513')

        ec.model = self.model

        self.measurer = f_me.ModelMeasurer(ec.Ve, self.model, ec.btrue) # сделали измерителя
        self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану

        #self.planner = f_p.DOptimalPlanner(ec, 'cache', verbose=True)

        self.planner = f_p.UniformPlanner(ec)

        #формирование цепочки
        estimator = f_e.NGEstimator()
        estimator.init_parameters(ec, self.model)

        ce  = f_e.ConsoleEstimatorDecorator(estimator)

        self.estimator = f_e.GraphPackEstimatorDecorator(ce,'../Cases/resfiles', 'NPN_KT315')

     # def proceed(self, nocacheplan = False):
     #    # Вот здесь бы и фабричный метод пользовать, потому что получается, что этот методв
     #    # требует некоторых атрибутов. Явно видно, какие методы надо переопределить,
     #    # а вот насчёт атрибутов не всё так прозрачно.
     #        plan = self.planner.give_us_a_plan(nocacheplan)
     #        self.measdata = self.plan_measurer.measure_acc_to_plan(plan)
     #
     #
     #        [print (m['x'], m['y']) for m in self.measdata]
     #
     #        exit(0)
     #
     #        return self.estimator.estimate(self.measdata, self.options)




def test():
    #dm = TransistorMainScript()
    dm = DiodeMainScript()


    #http://habrahabr.ru/post/157537/ - как накрутить производительность с помощью ctypes


    with f_sf.Profiler() as p:
        #print (dm.proceed(nocacheplan=False))
        dm.proceed(nocacheplan=False)


if __name__ == '__main__':
    test()

