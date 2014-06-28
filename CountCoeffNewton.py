__author__ = 'reiner'

import math

import numpy as np

import derivations as der


"""
В файле реализован метод подгонки коэффициентов методом Ньютона-Гаусса

"""



def readFile(filename):
    """
    читывает файл в массив списков (или векторов)
    """
    res=list()
    try:
        with open (filename, "r") as file:
            #print (file.readline())
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

    dic = {k:v for k,v in zip(totalstrlist, totallist)}   #with Python 3.x, goes for dict comprehensions   http://legacy.python.org/dev/peps/pep-0274/

    try:
        res=list()
        for outvar in outvarstrlist:
            res.append(eval(funcstrdict[outvar], dic))
    except BaseException as e:
       print ("error in countfunctvect:", e)
       print ("totalstrlist", totalstrlist)
       print ("invarlist", invarlist)
       print ("coefflist", coefflist)



    return res







def grandCountGN(funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, vrslst, NSIG=3):
    """
    funcstrdict - словарь строкового представления функций
    outvarstrlist -  список выходных переменных (y1, y2)
    invarstrlist - список входных переменных
    coeffstrlist - список коэффициентов (r1, r2, r3)

    #filename - имя файла с результатами эксперимента

    vrslst - вывод функции generate, данные эксперимента, в виде списка словарей
    NSIG=5 - количество значащих (точность)


    """
    log=""#строка, куда пишутся всякие сообщения

    if vrslst==None:
        print ("grandCountGN Error: cannot read file")
        return None
    #надо произвести два списка: список векторов Xs, и Ys из входного
    Xs=list()
    Ys=list()
#здесь можно ещё поиграть с лямбдами, чтоб полностью отказаться от итеративных  процессов
    for line in vrslst:
        la=lambda x: line[x]
        Xs.append(np.array (list (map (la, invarstrlist))))
        Ys.append(np.array (list (map (la, outvarstrlist))))



    #k=np.ones(len(coeffstrlist)) #начальное приближение вектора коэффициентов

    k=np.array((range (1, len(coeffstrlist)+1  ))   )

    k=np.array( [10, 20, 340])



    prevk=k #предыдущее значение вектора коэфф
    convergence=0
    numIterations=1



    A=np.zeros ((len(coeffstrlist), len(coeffstrlist)))
    b=np.zeros((len(coeffstrlist), 1))

    Sk=0
    Skmu=0
    N=len(Xs)  #размер выборки

    ind=0
    func=lambda x,k: np.array(countfunctvect (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, x.tolist(), k.tolist()))
    for xx in Xs:
        dif=Ys[ind]-np.array(func(xx,k))
        Sk+= np.dot(dif.T, dif)
        ind+=1



    Skpriv=0
    mu=1

    condition = True
    fstruct = lambda x,k: der.Jakobeand (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, x.tolist(), k.tolist())

    Tv=lambda x: (np.asmatrix(x)).T


    while condition:
        Skpriv=Sk
        prevk=k
        Sk=0
        A=np.zeros_like(A)
        b=np.zeros_like(b)




        for i in range (0, len(Xs)):
            fstructval=fstruct(Xs[i], k)
            A+=np.dot (fstructval.T, fstructval)
            ydif=Ys[i]-func(Xs[i],k)
            b+=np.dot (fstructval.T, Tv(ydif))   #транспонирование введено для согласования, не коррелирует с формулами


        deltak=np.linalg.solve(A,b)

        mu=2

        cond2=True
        it=0
        while (cond2):
            Skmu=0
            mu/=2
            for i in range (0, len (Xs)):

                vvv=Ys[i]-func(Xs[i], mu*deltak.T[0] + k)
                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)


                Skmu+=np.dot(vvv.T, vvv)


            it+=1
            if (it>1000):
                break
            cond2=Skmu>Skpriv

#        k+=mu*deltak
        k+=mu*deltak.T[0]
                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)




        Sk=Skmu


        numIterations+=1
        convergence=0

        for i in range (0, len (coeffstrlist)):
            convergence+=math.fabs(deltak[i]/prevk[i])
        convergence/=len(coeffstrlist)


        log+="Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n"



        print ("Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n")


        if (numIterations>100): #для ради безопасности поставим ограничитель на число итераций
            break
        condition = convergence>math.pow(10, -1*NSIG)


    #print (log)


    #пытаемся проанализировать результат: выводим средний остаток по функции с текущим K
    #по сути тестареа
    testdiff=0

    for i in range (0, len(Xs)):
        testdiff+=math.fabs(func(Xs[i], k)[1] - Ys[i][1])
    testdiff/=len(Xs)


    print ("testdiff: ", testdiff)


    return k, Sk, numIterations, testdiff



   # nvars = len(expdata[0])-len(funcstrstrlist) #количество входных переменных

    #print (nvars)




#test area


import sequence_generation as sg

funcstrdict= {"y1":"u1* (r2+r3)/(r1+r2+r3)", "y2":"u1* r3/(r1+r2+r3)"}
#funcstrdict= {"y1":"u1*r23/(r1+r23)", "y2":"u1*r3/(r1+r23)"}







xvectorlistsdict = {"u1":range(1,10),  "r1":[20], "r2":[30], "r3":[400]}

vrslst=sg.generate (funcstrdict, xvectorlistsdict, None, Vx=None, nvolx=None, yvectordispsdict=None, nvoly=1, outfilename="t.txt", listoutvars=["y1", "y2", "u1"] )

import pickle
#pickle.dump(vrslst, open("vrslt.f", "wb"))

vrslst=pickle.load(open("vrslt.f", "rb"))

#сюда впилить чтение файла

grandCountGN(funcstrdict,["u1"] , ["y1", "y2"],["r1", "r2", "r3"] , vrslst, NSIG=5)
#grandCountGN(funcstrdict,["u1"] , ["y1", "y2"],["r1", "r23"] , vrslst, NSIG=5)





