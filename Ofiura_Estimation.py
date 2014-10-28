__author__ = 'reiner'

import math
import copy

import numpy as np

"""
def grandCountGN_UltraX1 (funcf, jacf,  measdata:list, binit:list, c, NSIG=3):
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param c словарь дополнительных постоянных
    :param NSIG=3 точность (кол-во знаков после запятой)
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения


"""


def grandCountGN_UltraX1 (funcf, jacf,  measdata:list, binit:list, c, NSIG=3):
    """
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param c словарь дополнительных постоянных
    :param NSIG=3 точность (кол-во знаков после запятой)
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
    """

    b=binit
    log=""
    numiter=0
    condition=True
    while (condition):
        m=len(b) #число коэффициентов
        G=np.zeros((m,m))
        B5=None
        bpriv=copy.copy(b)
        Sk=0
        for point in measdata:
            jac=jacf(point['x'],b,c,point['y'])

            print (G.shape, jac.T.shape, jac.shape)
            G+=np.dot(jac.T,jac)

            dif=point['y']-funcf(point['x'],b,c)
            if B5==None:
                #B5=np.dot(jac,dif)
                B5=np.dot(dif, jac)

                print (dif.shape, jac.shape)

            else:
                B5+=np.dot(dif,jac)
            Sk+=np.dot(dif.T,dif)

        print (G.shape, '\n\n\n', B5.shape)

        deltab=np.dot(np.linalg.inv(G), B5)

        print ("Sk:",Sk)

        #mu counting
        mu=4
        cond2=True
        it=0
        while (cond2):
            Skmu=0
            mu/=2
            for point in measdata:
                dif=point['y']-funcf(point['x'],b-deltab*mu,c)
                Skmu+=np.dot(dif.T, dif)
            it+=1
            if (it>100):
                log+="Mu counting: break due to max number of iteration exceed"
                break
            cond2=Skmu>Sk

        b=b-deltab*mu

        print ("Iteration {0} mu={1} delta={2} deltamu={3} resb={4}".format(numiter, mu, deltab, deltab*mu, b))

        numiter+=1

        condition=False
        for i in range (len(b)):
            if math.fabs ((b[i]-bpriv[i])/bpriv[i])>math.pow(10,-1*NSIG):
                condition=True

        if numiter>100: #max number of iterations
            log+="GKNUX1: Break due to max number of iteration exceed"
            break


    return b, numiter, log
