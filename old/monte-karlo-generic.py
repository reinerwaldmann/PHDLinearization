__author__ = 'vasilev_is'

import numpy as np
import matplotlib.pyplot as plt
import random


def evaldfunc (funcstr, arginitseq1):
    """
    вычисляет функцию,
    arginitseq1 - locals
    funcstr - function line


    """
    return  eval(funcstr, arginitseq1)



#генерирует список случайных чисел м=0 д=1
def generlist (n):
    r=list()
    for x in range (0,n):
        r.append(random.gauss(0,1))
    return r



def generrandvals (M, cov_ksi, nvol=1000):
    """
    генерирует кортеж случайных последовательностей, коррелированных по матрице и с матожиданием во вх. параметрах
    """
    M_ksi=np.array(M)
    #проверка корректности:

    if (np.linalg.det (cov_ksi) <0 ):
        print ("Определитель матрицы меньше  0")
        return None


    n=len(M_ksi)  #размерность вектора
    U=list() #генерация случайного вектора U


    for x in range (0,nvol):
        U.append(generlist(n))

    A=np.linalg.cholesky(cov_ksi)

    ksi=list()
    for i in range (0, nvol):
        ksi.append( np.dot(A, np.array(U[i]).T)+M_ksi.T)

    #ksi - список значений случайных векторов, значения векторов коррелированы (значение - вектор)

    #эта часть разбрасывает значения аккуратно по спискам, так, что получается список векторов (вектор - значение)

    res=list()
    for i in range (0, n):
        res.append(list())

    for ni in range (0,nvol):
        for i in range (0, n):
            res[i].append (ksi [ni][i])

    XX=np.array(res)
    COVTEST=np.cov(XX, bias=-1)

    MTEST=list(map (np.mean, res))


    # tseq1=list()
    # tseq2=list()
    # tseq3=list()
    # for ni in range (0,nvol):
    #     tseq1.append(ksi[ni][0])
    #     tseq2.append(ksi[ni][1])
    #     tseq3.append(ksi[ni][2])
    # XX=np.array([tseq1, tseq2, tseq3])
    # COVTEST=np.cov(XX, bias=-1)

    #print ("Mean Values:")
    #print (np.mean(tseq1), np.mean(tseq2), np.mean(tseq3)  )
   # print ("Cov matrix:")
    #print (COVTEST)
    #plt.hist(tseq1, 100, label="tseq1")
    #plt.xlabel('Smarts')
    #plt.ylabel('Probability')
    #plt.title("TSEQ1")
    #plt.show()
    return res, COVTEST, MTEST


#считает дисперсию и матожидание выхода, которые получаются, ежели применить сгенерированные последовательности к функции



def countDispMonteKarlo (fun_seq, argseq, arginitseq1, V, nvol=1000):

    M=list()
    for i in range (0, len(argseq)):
        M.append(arginitseq1[argseq[i]])

    #теперь в M матожидания располагаются так же, как расположены переменные argseq


    #получаем список изслучайных векторов, которые коррелированы, как нам надо
    randvalsvect=generrandvals (M, V, nvol)


    #печатаем всё, чтоб было ясно, что всё ОК
    # print (argseq)
    # print (randvalsvect[2])
    # print (randvalsvect[1])
    # print (randvalsvect[0])



    #список значений функции при случайных данных
    listfunresvals=dict()

    varfunc=dict()
    meanfunc=dict()

    for fun in fun_seq: #для каждой функции из списка функций
        listfunresvals[fun]=list()   #инициализируем список результатов функции
        for i in range (0, nvol):  #для каждого значения из случайных списков
            #формируем значения переменных: есть список переменых с разбросом, есть три списка их
            # случайных значений, есть словарь, который надо обновить значениями переменных с разбросом
            arginitseq1here=dict()
            j=0
            for k in argseq:
                arginitseq1here[k]=randvalsvect[0][j][i]


                j += 1
                #ainitsseq=arginitseq1
                #ainitsseq.update(arginitseq1here)  #добавляем значения переменных без разброса

            arginitseq1here["u1"]= 100
            listfunresvals[fun].append (eval (fun, arginitseq1here))




            #добавляем новоее значение  функции

        meanfunc[fun]=np.mean(listfunresvals[fun])
        varfunc[fun]=np.var(listfunresvals[fun])

    return varfunc, meanfunc




def countDispMonteKarloWrapper (fun_seq, argseq, arginitseq1, V, nvol=1000):
    """
    Контролирует правильность аргументов
    """
    return countDispMonteKarlo (fun_seq, argseq, arginitseq1, V, nvol)



#testing area

M=np.array([20,30,400])
funcseq= ["u1* (r2+r3)/(r1+r2+r3)", "u1* r3/(r1+r2+r3)"]



argseq=["r1", "r2", "r3" ]
arginitseq1={"u1":100, "r1":M[0], "r2":M[1], "r3":M[2]}

V=np.array       ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])

V1=np.array      ( [[4, 0, 0],
                    [0, 9, 0],
                    [0, 0, 16]])


res=countDispMonteKarloWrapper(funcseq, argseq, arginitseq1, V, 1000)
print (res[0])
print (res[1])

