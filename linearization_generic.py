__author__ = 'reiner'

#Все формулы получены по статье про транзисторы и
# Теория вероятностей.  Вентцель Е.С.
#4-е изд., стереотип. - М.: Наука, Физматгиз, 1969 - 576 с.

import numpy as np
import sympy as smp

def derivsym (funcstr, argseq, n=1):
    """
    Производит символьное дифференцирование входящей функции по переменным из последовательности argseq
    возвращает словарь переменная - строковое представление функции
    funcstr - строка, являющая собою функцию
    argseq - последовательность из аргументов, заданных как строки
    n - номер производной

    """
    res=dict()
    for arg in argseq:
        res[arg]=str(smp.diff(funcstr,arg, n))
    return res

def evalderivsymv (funcstr, argseq, arginitseq, n=1):
    """
    Вычисляет значения производных функции по всем аргументам, заданным в argseq
    arginitseq - значения переменных (всех)
    funcstr - строка, являющая собою функцию
    argseq - последовательность из аргументов, заданных как строки
    возвращает словарь: переменная - значение производной при приложении словаря значений переменных arginitseq
    """
    res=dict()
    derivdict=derivsym (funcstr, argseq, n)
    for arg in argseq:
        res[arg] = eval(derivdict[arg], arginitseq)
    return res


def deriv_func_seq (fun_seq, argseq, arginitseq, n=1):
    """
    Производит символьное дифференцирование последовательности функций. Возвращает значения производных
    На вход получает список из строковых значений функций,
    на выходе возвращает словарь строковое функции - значение производной
    при приложении arginitseq, как словаря значений переменных.
    n - номер производной
    """
    res=dict()
    for fun in fun_seq:
        res[fun] = evalderivsymv(fun,argseq, arginitseq, n)
    return res



def countDispLinearization (fun_seq, argseq, arginitseq, V):
    """
    Cчитает дисперсию по матожиданию и ковариационной матрице методом линеаризации

    fun_seq - последовательность функций, заданных строками
    argseq - список аргументов в той последовательности, в каковой они расположены в ковариационной матрице
    arginitseq - словарь переменная - значения. Здесь должны быть поименованы все переменные, даже не имеющие разброса
    #V - корреляционная матрица
    """


    derivdict = deriv_func_seq(fun_seq, argseq, arginitseq, 1)

    res=dict()

    #для каждой функции из списка функций
    for fun in fun_seq:
        derivfunc = derivdict[fun]  #таким образом мы получили словарь значений производных
        first_member=0
        for i in range (0, len(argseq)):
            first_member+=(derivfunc[argseq[i]]**2) * V[i,i]

        second_member=0
        for i in range (0, len(argseq)):
            for j in range (0, len(argseq)):
                if (i<j):
                    second_member += 2*derivfunc[argseq[i]]*derivfunc[argseq[j]]*V[i,j]

        res[fun]=first_member+second_member

    return res




def countDispLinearizationWrapper (fun_seq, argseq, arginitseq, V):
    """
    Контролирует правильность аргументов
    """
    return countDispLinearization (fun_seq, argseq, arginitseq, V)


#testing area

M=np.array([20,30,400])
funcseq= ["u1* (r2+r3)/(r1+r2+r3)", "u1* r3/(r1+r2+r3)"]



argseq=["r1", "r2", "r3" ]
arginitseq={"u1":100, "r1":M[0], "r2":M[1], "r3":M[2]}

V=np.array       ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])

V1=np.array      ( [[4, 0, 0],
                    [0, 9, 0],
                    [0, 0, 16]])


print (countDispLinearizationWrapper(funcseq, argseq, arginitseq, V))


#import resistors
#resistors.test2(M,V,1000)


