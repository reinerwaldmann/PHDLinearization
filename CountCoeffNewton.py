__author__ = 'reiner'

import numpy as np

def readFile(filename):
    """
    читывает файл в массив списков (или векторов)
    """
    res=list()
    try:
        with open (filename, "r") as file:
            print (file.readline())
            for line in file:
                res.append (list(map (lambda x: float(x), line.split("\t")[:-1])))

        return res
    except BaseException:
        return None


def countfunctvect (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, invarlist, coefflist):
    """
    Считает значение функции при приложении вектора значений входов
    funcstrdict строкое пр функций

    invarstrlist список названий входных перем
    outvarstrlist список названий вых пер
    coeffstrlist список названий коэфф

    invarlist список вх перем
    coefflist список коэфф
    """
    totalstrlist = invarstrlist+coeffstrlist
    totallist = invarlist+coefflist
    dic = {k:v for k,v in zip(totalstrlist, totalstrlist)}   #with Python 3.x, goes for dict comprehensions   http://legacy.python.org/dev/peps/pep-0274/

    res=list()
    for outvar in outvarstrlist:
        res.append(eval(funcstrdict[outvar],locals=dic))


    return res







def grandCountGN(funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, vrslst, NSIG=5):
    """
    funcstrdict - словарь строкового представления функций
    outvarstrlist -  список выходных переменных (y1, y2)
    invarstrlist - список входных переменных
    coeffstrlist - список коэффициентов (r1, r2, r3)

    #filename - имя файла с результатами эксперимента

    vrslst - вывод функции generate, данные эксперимента, разложенные в словарь списков

    NSIG=5 - количество значащих (точность)


    """
    log=""#строка, куда пишутся всякие сообщения

    if vrslst==None:
        print ("grandCountGN Error: cannot read file")
        return None

    k=np.ones(len(coeffstrlist)) #начальное приближение вектора коэффициентов
    prevk=k #предыдущее значение вектора коэфф
    convergence=0
    numIterations=1

    A=np.zeros ((len(coeffstrlist), len(coeffstrlist)))
    b=np.zeros((len(coeffstrlist, 1)))

    Sk=0
    Skmu=0

    #N=len(vrslst.values()[0])  #размер выборки







   # nvars = len(expdata[0])-len(funcstrstrlist) #количество входных переменных

    #print (nvars)










    pass




#test area


import sequence_generation as sg

funcstrdict= {"y1":"u1* (r2+r3)/(r1+r2+r3)", "y2":"u1* r3/(r1+r2+r3)"}

xvectorlistsdict = {"u1":range(1,10),  "r1":[20], "r2":[30], "r3":[400]}

vrslst=sg.generate (funcstrdict, xvectorlistsdict, None, Vx=None, nvolx=None, yvectordispsdict=None, nvoly=1, outfilename="t.txt", listoutvars=["y1", "y2", "u1"] )






