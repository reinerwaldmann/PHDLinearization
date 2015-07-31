
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

    def proceed(self):
            #прокатит и за абстрактный метод

            plan = self.planner.give_us_a_plan()
            self.measdata = self.plan_measurer.measure_acc_to_plan(plan)
            # естиматор инициализируется как-то уж слишком через задницу - в настоящее время
            # идёт позним связыванием - типа поле объекта лишь ссылка.
            self.estimator.estimate(self.measdata, f_e.Options())




class DiodeMainScript(AbstractMainScript):
    def __init__(self):

        btrue = [7.69e-8, 1.45 ,.0422] #номинальные значения диода D1N4001 с сайта, вроде официальной модели производителя
        binit = btrue
        Ve=np.array([[1.9e-5]])
        bstart=np.array(self.btrue)-np.array(self.btrue)*0.2
        bend=np.array(self.btrue)+np.array(self.btrue)*0.2
        xstart=[0.001]
        xend=[1]
        N=20
        ec =f_sf.EstimationContext(bstart, bend, btrue, binit, xstart, xend, Ve, N)



        self.model = f_m.SimpleDiodeModel ('Diode_1N') # сделали модель

        self.measurer = f_me.ModelMeasurer(self.Ve, self.model, self.btrue) # сделали измерителя

        self.plan_measurer = f_me.PlanMeasurer(self.measurer) # сделали измеритель по плану

        #self.planner = f_p.DOptimalPlanner(self.N, self.xstart, self.xend, self.bstart, self.bend, self.model, verbose=1)
        self.planner = f_p.UniformPlanner(self.N, self.xstart, self.xend)


        self.measdata=None

        #формирование цепочки
        estimator = f_e.NGEstimator()
        estimator.init_parameters(self.bstart, self.bend, self.binit, self.xstart, self.xend, self.model, self.measdata, self.Ve)

        self.estimator  = f_e.ConsoleEstimatorDecorator(estimator)




def test():
    dm = DiodeMainScript()
    dm.proceed()


test()