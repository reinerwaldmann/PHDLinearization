__author__ = 'reiner'

import hashlib

import sympy as smp
import numpy as np

"""
В файле реализованы функции символьного диффенцирования векторных функций с оптимизацией

"""

derivdict=dict()

def makefunnyhash (funcstr, n):
    """
    Делает имя файла из функции и номера производной.
    Требует наличия доступной папки funcs для работы скрипта, иначе производит ошибки
    """

    return  "funcs/"+hashlib.sha224(funcstr.encode()).hexdigest()[0:30]+str(n)

    #return "funcs/"+funcstr+"_"+str(n)




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


    funnyhash=makefunnyhash (funcstr+argseq.__str__(), n)

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
            filef.write(arg+"\t"+res[arg]+"\n")


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

    global derivdict

    funnyhash = makefunnyhash (funcstr+argseq.__str__(), n)  #определить требуемое имя файла
    if not (funnyhash in derivdict):
        try:
            with open (funnyhash, 'rt') as file:  #попытаться открыть файл, где записаны производные функции по переменной. Одна функция - один файл!
                thisfundict=dict()
                for line in file:
                    l=line.split("\t")
                    thisfundict[l[0]]=l[1]  #считать оттуда символные представления производных
                derivdict[funnyhash] = thisfundict

        except BaseException: #если не вышло с файлом
            derivdict[funnyhash]=derivsym (funcstr, argseq, n) #то получить символьные значения производных (оно и файл запишет)



    for arg in argseq:

        res[arg] = eval(derivdict[funnyhash][arg], arginitseq) #вычислить значения производных
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
        res[fun] = evalderivsymv(fun, argseq, arginitseq, n)


    return res

def Jakobean (fun_seq, argseq, arginitseq):
    """
    Возвращает якобиан векторной функции, по сути матричная форма deriv_func_seq

    fun_seq - последовательность строковых записей функций в должном порядке
    argseq - последовательность аргументов с разбросом в должном порядке
    arginitseq - значения всех аргументов

    Якобиан
    это функция, возвращающая np.array (то есть, матрицу по сути),
    являющая собой производные всех функций по всем переменным,
    df2/dfx1  df2/dx2
    df1/dfx1  df1/dx2
    """

    res = np.zeros ( (len(fun_seq), len(argseq))  )


    fndrv=deriv_func_seq (fun_seq, argseq, arginitseq, n=1)

    #print (fndrv)

    for i in range (0, len(fun_seq)):
        for j in range (0, len(argseq)):
            res[i][j]=fndrv[fun_seq[i]][argseq[j]]

#            print (fun_seq[i], argseq[j])



    return res


def Jakobeand (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, invarlist, coefflist):

    flst=list(map(lambda x:funcstrdict[x],outvarstrlist))

    totalstrlist=invarstrlist+coeffstrlist
    totallist=invarlist+coefflist
    dic = {k:v for k,v in zip(totalstrlist, totallist)}   #with Python 3.x, goes for dict comprehensions   http://legacy.python.org/dev/peps/pep-0274/



    return Jakobean(flst, coeffstrlist, dic)









#test area
#тестируем на время: в первый раз 0.02 сек, далее 0.0 сек, когда файл уже создан, т. к. считывает сначала из файла, потом вообще из кеша.

# from time import clock
#
# for i in range (0,10):
#     start1 = clock()
#     test=evalderivsymv ("2*x+1+f*x", ['x'], {'x':1 ,'f':100} )
#     end1 = clock()
#     print("Result (iterativ): ", test, "\nDie Funktion lief %1.10f Sekunden" % (end1 - start1))

funcstrdict= {"y1":"u1* (r2+r3)/(r1+r2+r3)", "y2":"u1* r3/(r1+r2+r3)"}
xvectorlistsdict = {"u1":100,  "r1":20, "r2":30, "r3":400}
spreadvarslist  = ["r1", "r2", "r3"]







#derivsym ("2*x+1", ('x'), 1)