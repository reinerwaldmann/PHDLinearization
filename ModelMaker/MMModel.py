__author__ = 'vasilev_is'


import sympy as smp
import numpy as np


class MMModel:
    """
    Описывает модель.
    Это однооткликовая версия, потому многие переменные (как  _modelstring идут как строка, а должны бы как список)
    """

    def __init__(self):
        self.assigns=[]  #список назначений рода x1=x[0]
        self.modelparts=[]  #список суммируемых членов модели

    #    self._modelstring=''
     #   self._jacstring=''

    def showModel(self):
        strmodel = 'y0='+'+'.join(self.modelparts) #возвращает строку модели
        strassigns = '\n'.join(self.assigns)
        evalstr = strassigns+'\n'+strmodel

        print (evalstr)






    def solver (self, x,b,c=None):
        strmodel = 'y0='+'+'.join(self.modelparts) #возвращает строку модели
        strassigns = '\n'.join(self.assigns)
        evalstr = strassigns+'\n'+strmodel

        ns={}

        y0=None

        exec(evalstr, locals(), ns)
        # немного об exec в данном случае. Она принимает на вход, кроме строки кода, ещё массивы globals и locals
        # цимес в том, что она их фиксирует в момент выполнения, то есть изменить локальные и глобальные
        # переменные exec не может чтобы получить из неё параметры, нужно один из этих массивов изменить на
        # пользовательский. К примеру, здесь вместо глобальных подсунуты локальные, а вместо локальных - пустой словарь.
        # В нём, по итогам, и окажутся все глобальные переменные. Есть, однако, некоторые сомнения в том, а
        # хорошо ли так делать. Возможно, лучше пойти некоторым другим путём.


        return [ns['y0']]


    def jacf (self,  x,b,c,y):
        """
        Формирует и возвращает якобиан по коэффициентам.
        :param x:
        :param b:
        :param c:
        :param y:
        :return:
        """
        strmodel = '+'.join(self.modelparts) #возвращает строку модели, её будем дифференцировать
        listbletters = ['b{0}'.format(i)  for i in range(len(b))] #получаем список вида b0, b1 и так далее
        bdiffstrlst = [str(smp.diff(strmodel, b1)) for b1 in listbletters ] #получили производную
        strassigns = '\n'.join(self.assigns)
        diflist=[]

        for dif in bdiffstrlst:
            evalstr = strassigns+'\n'+'d0='+dif



            ns={}
            exec(evalstr, locals(), ns)
            diflist.append(ns['d0'])

        return np.matrix([diflist])


    def makeLinearRegression(self, nvar):
        """
    Функция уровнем повыше: создать уравнение линейной регрессии
        :param nvar: число переменных регрессии
        :return:
        """
        self.assigns = ['x{0}=x[{0}]'.format(i) for i in range (nvar)]
        self.assigns+= ['b{0}=b[{0}]'.format(i) for i in range (nvar+1)]

        self.modelparts = ['b0']
        self.modelparts += ['b{0}*x{1}'.format(i+1, i) for i in range(nvar) ]
        #http://habrahabr.ru/post/63539/
        # это пост на хабре про сложение списков. Он и говорит, что самый быстрый способ -
        # имено такой - extend или +=

        # http://pythonworld.ru/tipy-dannyx-v-python/vse-o-funkciyax-i-ix-argumentax.html
        # статья о функциях с произвольным числом аргументов



def test():
#    assigns=['x0=x[0]', 'x1=x[1]', 'x2=x[2]'] +  ['b0=b[0]', 'b1=b[1]', 'b2=b[2]', 'b3=b[3]']
 #   modelparts=['b0', 'b1*x0', 'b2*x1', 'b3*x2']
    x=[1,2,3]
    b=[1,2,3,4]
    #
    # model = MMModel()
    #
    # model.assigns = assigns
    # model.modelparts = modelparts
    #
    # print (model.jacf(x,b,None, None))
    model = MMModel()
    model.makeLinearRegression(3)
    model.showModel()
    print (model.solver(x,b))
    print (model.jacf(x,b,None,None))

#test()

