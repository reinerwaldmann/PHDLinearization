__author__ = 'vasilev_is'

import random
import math
import numpy as np


class AbstractMeasurer:
    "Defining a measurer abstraction"
    #_xstart=[]
    #_xend=[]

    def measure (self, x):
        raise NotImplementedError("Please Implement this method")

    def getCovMatrix(self):
        raise NotImplementedError("Please Implement this method")


class ModelMeasurer(AbstractMeasurer):

    def __init__(self, Ve, model, b):
        self.Ve = Ve
        self.model = model
        self.b = b

    def measure(self, x):
        y=self.model.funcf(x,self.b)
        if y is None: #если функция вернула чушь, то в measdata её не записывать!
            return None
        #Внесём возмущения:
        if self.Ve is not None and np.linalg.det(self.Ve)>10e-15:
            ydisps=np.diag(self.Ve)


            for k in range(len(y)):
                if (y[k] >=  0):
                    y[k]=math.exp(random.normalvariate(math.log(y[k]), math.sqrt(ydisps[k])))
                elif (y[k] < 0):
                    y[k]=math.exp(random.normalvariate(-1*math.log(math.fabs(y[k])), math.sqrt(ydisps[k])))


        return y

    def getCovMatrix(self):
        return self.Ve



class PlanMeasurer:
    "Performing measurements according to a plan"
    def __init__(self, measurer):
        self.measurer = measurer

    def measure_acc_to_plan (self, plan):
        res=[]
        for point in plan:
            res.append({'x':point, 'y':self.measurer.measure(point)})
        return res

