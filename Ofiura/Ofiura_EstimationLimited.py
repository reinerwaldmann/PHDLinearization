__author__ = 'vasilev_is'

import math
import copy
import platform

import numpy as np


"""
ПОстроена на использовании библиотеки numpy
"""




def makeSkInit (funcf, measdata, binit, c):
    Sk=None #или mpm.mpf(0.0)
    for point in measdata:
        fxbc=funcf(point['x'],binit,c)
        if fxbc is None:
            print ("ERROR: makeSkInit: funcf returned None")
            return None
        dif=point['y']-fxbc
        if Sk is None:
            Sk=np.dot(dif.T,dif)
        else:
            Sk+=np.dot(dif.T,dif)
    return Sk

def islinux ():
     if 'Windows' in platform.system():
         return False
     return True

def makeAinit (bstart, bend, Skbinit, binit, isBinitGood=True):
    """
    расчёт ведётся по правилу "binit есть хорошее предположение в рамках области"
    :param bstart:
    :param bend:
    :param Skbinit:
    :param binit:
    :return:
    """
    A=np.zeros ( (len(binit), 2 ))  #,  dtype=np.float96 if islinux() else np.float64

    for i in range (len(binit)):
        #A[i][0]=A[i][1]=.001*(bend[i]-bstart[i])*Skbinit  #универсально

        if isBinitGood:
            A[i][0] = 0.001*(binit[i]-bstart[i])*Skbinit
            A[i][1] = 0.001*(bend[i]-binit[i])*Skbinit
        else:
            A[i][0] = A[i][1] = 0.001*(bend[i]-bstart[i])*Skbinit  #универсально
    return A

def countSklims(A,b,bstart, bend):
    partone=parttwo=.0
    for i in range (len(b)):
        partone+=A[i][0]/(b[i]-bstart[i]) #UNSAFE: если вдруг, ну вдруг b=bstart или bend, то будет <strike> треш </strike> креш с делением на ноль.
        parttwo+=A[i][1]/(bend[i]-b[i])
    return partone+parttwo

def countN (A, b, bstart, bend):

    N=np.zeros ((len(b),len(b)))
    for j in range (len(b)):
        partone=2*A[j][0]/(b[j]-bstart[j])**3
        parttwo=2*A[j][1]/(bend[j]-b[j])**3
        N[j][j]+=parttwo+partone #так как матрица нулевая
    return N

#TODO требуется добавить функцию с оценкой качества, вопрос, сюда ли это добавить или же сделать враппер-функцию оценки качества в qualitat файле, и её надевать на разные gknux


def  grandCountGN_UltraX1_Limited_wrapper (funcf, jacf,  measdata:list, binit:list, bstart:list, bend:list, c, NSIG=50, NSIGGENERAL=50, implicit=False, verbose=False, verbose_wrapper=False, isBinitGood=True):
    """
    Обёртка для grandCountGN_UltraX1_Limited для реализации общего алгоритма
    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param bstart:list нижняя граница b
    :param bend:list верхняя граница b
    :param c словарь дополнительных постоянных
    :param A матрица коэффициентов a
    :param NSIG=3 точность (кол-во знаков после запятой)
    :param implicit True если функция - неявная, иначе false
    :param verbose Если True, то подробно принтить результаты итераций
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
    """

    maxiter=100
    b,bpriv=binit,binit
    gknux=None

    Skinit = makeSkInit(funcf,measdata,binit, c)
    A=makeAinit(bstart, bend,Skinit,binit, isBinitGood)

    log=''

    if verbose_wrapper:
        print ('==grandCountGN_UltraX1_Limited_wrapper is launched==\n\n')

    for numiter in range (maxiter):
        bpriv=copy.copy(b)
        gknux=grandCountGN_UltraX1_Limited (funcf, jacf,  measdata, b, bstart, bend, c, A, NSIG, implicit, verbose) #посчитали b

        if verbose_wrapper:
            print ('Iteration \n',numiter,'\n' ,gknux)
        b=gknux[0]

        if not gknux[2]=='':
            log+=gknux[2]+"\n"
        for j in range (len(binit)): #уменьшили в два раза
            A[j][0]*=0.5
            A[j][1]*=0.5

        condition=False
        for i in range (len(b)):
            if math.fabs ((b[i]-bpriv[i])/bpriv[i]) > math.pow(10,-1*NSIGGENERAL):
                condition=True
        if not condition:
            break

        #мол если хоть один компонент вектора b значимо изменился, тогда продолжать. Иначе программа дойдёт до break и цикл прекратится

#  return b, numiter, log, Sklist, Sk

    print ('grandCountGN_UltraX1_Limited_wrapper iterations number:', numiter)
    return gknux[0], gknux[1], log, gknux[3], gknux[4]




def grandCountGN_UltraX1_Limited (funcf, jacf,  measdata:list, binit, bstart, bend, c, A, NSIG=3, implicit=False, verbose=False):
    """
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом с ограничениями, заданными диапазоном
    В стандартный поток вывода выводит отладочную информацию по каждой итерации

    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param bstart:list нижняя граница b
    :param bend:list верхняя граница b
    :param c словарь дополнительных постоянных
    :param A матрица коэффициентов a
    :param NSIG=3 точность (кол-во знаков после запятой)
    :param implicit True если функция - неявная, иначе false
    :param verbose Если True, то подробно принтить результаты итераций
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
    """
    #sign - если  1, то b=b+deltab*mu, иначе b=b-deltab*mu. При неявной функции надо ставить sign=0
    sign=0 if implicit else 1

    Sklist=list()
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

            if jac is None:
                return None

            #print (G.shape, jac.T.shape, jac.shape)
            G+=np.dot(jac.T,jac)

            try:
                dif=np.array(point['y'])-np.array(funcf(point['x'],b,c))
            except BaseException as e:
                print('grandCountGN_UltraX1_limited: As funcf returned None, method  stops:', e)
                print (point, b, sep='\n')
                return None

            # print(dif, jac)
            # print('-----')
            #
            # print()
            if B5 is None:
                B5=np.dot(dif, jac)
            else:
                B5+=np.dot(dif,jac)
            Sk+=np.dot(dif.T,dif)


        N=countN(A,b,bstart, bend)
        Sklims=countSklims(A,b,bstart, bend)
       #G=G-N if implicit else G+N #добавляем градиент от штрафных функций
        G=G+N
        Sk+=Sklims #добавиляем объектную функцию от штрафных функций

        #print(np.linalg.inv(G), B5[:,0])
        #костыль для диодной задачи
        if hasattr(B5, 'A1'):
            B5=B5.A1
        try:
            deltab=np.dot(np.linalg.inv(G), B5)
        except BaseException as e:
            print('Error in G:', e)
            print('G=',G)
            return None


        #mu counting
        mu=4
        cond2=True
        it=0
        Skmu=0
        while (cond2):
            Skmu=countSklims(A,b,bstart, bend) #добавиляем объектную функцию от штрафных функций
            mu/=2
            for point in measdata:
                try:
                    dif=np.array(point['y'])-np.array(funcf(point['x'],b+deltab*mu,c)) if sign else np.array(point['y'])-np.array(funcf(point['x'],b-deltab*mu,c))
                except:
                    continue

                Skmu+=np.dot(dif.T, dif)

            it+=1


            if (it>100):
                log+="Mu counting: break due to max number of iteration exceed"
                break
            cond2=Skmu>Sk

        b=b+deltab*mu if sign else b-deltab*mu

        Sk=Skmu

        Sklist.append(Sk)
        if verbose:
            print ("Sk:",Sk)
            print ("Iteration {0} mu={1} delta={2} deltamu={3} resb={4}".format(numiter, mu, deltab, deltab*mu, b))

        numiter+=1

        condition=False
        for i in range (len(b)):
            if math.fabs ((b[i]-bpriv[i])/bpriv[i])>math.pow(10,-1*NSIG):
                condition=True


        if numiter>2000: #max number of iterations
            log+="GKNUX1: Break due to max number of iteration exceed"
            break


    return b, numiter, log, Sklist, Sk


#693175 усадьба сервис