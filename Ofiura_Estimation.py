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

            #print (G.shape, jac.T.shape, jac.shape)
            G+=np.dot(jac.T,jac)

            dif=np.array(point['y'])-np.array(funcf(point['x'],b,c))
            if B5==None:

                B5=np.dot(dif, jac)

                #print (dif.shape, jac.shape)

            else:
                B5+=np.dot(dif,jac)


            Sk+=np.dot(dif.T,dif)





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




                dif=np.array(point['y'])-np.array(funcf(point['x'],b-deltab*mu,c))
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


def grandCountGN_UltraX (funcf, jacf,  expdatalist:list, kinit:list, NSIG=3):
    """
    Подгоняет коэфф. методом Ньютона-Гаусса с изменяемой длиной шага (через mu)
    Parameters:
    funcf - функция, на вход которой подаются вектора x и b, а на выходе получается вектор y, притом без возмущений
    jacf - функция, на вход которой подаются вектора x и b, а на выходе получается якобиан функции
    expdatalist:list - экспериментальные данные
    kinit=None - начальное приближение коэффициентов
    NSIG=3 - число значащих цифр (точность подгонки коэффициентов)
    """
    log=""#строка, куда пишутся всякие сообщения

    if expdatalist==None:
        print ("grandCountGN_Ultra Error: cannot read exp data")
        return None
    #надо произвести два списка: список векторов Xs, и Ys из входного
    k=kinit
    M=len(k) # число оцениваемых коэффициентов
    prevk=k #предыдущее значение вектора коэфф
    convergence=0
    numIterations=1
    Sk=0
    Skpriv=0
    N=len(expdatalist)  #размер выборки
    condition = True
    while condition: #пока не пришли к конвергенции
        Skpriv=Sk
        prevk=k
        Sk=0
        PP=np.zeros ((M, M))
        PYY=np.zeros((M, 1))
        Tv=lambda x: (np.asmatrix(x)).T
        for i in range(0, N):


            PP+=np.dot(jacf(expdatalist[i]['x'], k, None).T, (jacf(expdatalist[i]['x'], k, None))  )

            dif = np.array(expdatalist[i]['y'])-np.array(funcf(expdatalist[i]['x'],k))
            Sk+= np.dot(dif.T, dif)
            PYY+=np.dot(jacf(expdatalist[i]['x'], k, None), np.transpose(np.asmatrix(dif)))
        deltak=np.dot(np.linalg.inv(PP), PYY)

        #применение mu
        mu=4
        cond2=True
        it=0
        while (cond2):
            Skmu=0
            mu/=2
            for i in range (0, len (expdatalist)):
                dif = np.array(expdatalist[i]['y'])-np.array(funcf(expdatalist[i]['x'], k+mu*deltak.T[0]))
                Skmu+=np.dot(dif.T, dif)
            it+=1
            if (it>100):
                print ("break")
                break
            cond2=Skmu>Sk

        k+=mu*deltak.T[0]

        Sk=Skmu



        numIterations+=1
        convergence=0

        for i in range (0, M):
            convergence+=math.fabs(deltak[i]/prevk[i])
        convergence/=M


        log+="Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n"

        print ("Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n")


        if (numIterations>100): #для ради безопасности поставим ограничитель на число итераций
            break
        condition = convergence>math.pow(10, -1*NSIG)


    #print (log)


    #пытаемся проанализировать результат: выводим средний остаток по функции с текущим K
    #по сути тестареа
    testdiff=0

    for i in range (0, len(expdatalist)):
        testdiff+=math.fabs(funcf(expdatalist[i]['x'], k)[1] - expdatalist[i]['y'][1]) #TODO проверить целесообразность [1]
    testdiff/=len(expdatalist)


    print ("testdiff: ", testdiff)


    #return k, Sk, numIterations, testdiff
    return k, numIterations, log