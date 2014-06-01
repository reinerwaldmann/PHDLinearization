__author__ = 'reiner'

import numpy as np
import sympy as smp
import hashlib


def derivsym (funcstr, argseq, n=1):
    """
    Производит символьное дифференцирование входящей функции по переменным из последовательности argseq
    возвращает словарь переменная - строковое представление функции
    funcstr - строка, являющая собою функцию
    argseq - последовательность из аргументов, заданных как строки
    n - номер производной
1
    """
    res=dict()
    for arg in argseq:
        res[arg]=str(smp.diff(funcstr,arg, n))

    funnyhash=hashlib.sha224(funcstr.encode()).hexdigest()[0:20]+str(n)+".py"

    try:
        open(funnyhash, 'r')

    except IOError:
        filef=open(funnyhash, 'wt')
    #нет никакой защиты от кривого открытия, шо пичалька


        # filef.write("def "+funnyhash+"(funcstr, arginitseq):")
        #
        # funcstring="\n\tres=dict()\n"
        # funcstring+="\tderivdict=["
        # for arg in argseq:
        #     funcstring+="'"+arg+"':'"+res[arg]+"',"
        # funcstring=funcstring[0:-1]+"]\n"
        #
        # funcstring += "\tfor arg in argseq: \n\t\tres[arg] = eval(derivdict[arg], arginitseq)\n\treturn res"
        # filef.write(funcstring)

        for arg in argseq:
            filef.write(arg+"\t"+res[arg])



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

    derivdict=dict()
    funnyhash=hashlib.sha224(funcstr.encode()).hexdigest()[0:20]+str(n)
    try:
        file=open (funnyhash, 'rt')
        for line in file:
            l=line.split("\t")
            derivdict[l[0]]=l[1]


    except BaseException:
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

def Jakobean (fun_seq, argseq, arginitseq):
    """
    Возвращает якобиан векторной функции, по сути матричная форма deriv_func_seq
    """
    pass



#TODO Jakobean
#TODO кеширование функции


#test area


print (evalderivsymv ("2*x+1+f*x", ['x'], {'x':1 ,'f':100} ))


#derivsym ("2*x+1", ('x'), 1)