

__author__ = 'vasilev_is'

"""
Назначение файла - тестирование возможности оценки параметров транзистора с помощью метода Ньютона-Гаусса
Функции - метод Ньютона-Гаусса в исходном виде, функция, якобиан, априорное планирование

Сперва возьмём обычную резистивную модель и попробуем её оценить.  С помощью обычного (равномерного) планирования, с помощью априорного плана
"""

import math

import numpy as np
from scipy import optimize


def grandCountGN_Ultra (funcf, jacf,  expdatalist:list, kinit:list, c, NSIG=3):
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
    A=np.zeros ((M, M))
    b=np.zeros((M, 1))
    Sk=0
    Skmu=0
    N=len(expdatalist)  #размер выборки

    func = funcf
    for i in range(0, len(expdatalist)):
        dif = np.array(expdatalist[i]['y'])-np.array(func(expdatalist[i]['x'],k,c))
        Sk+= np.dot(dif.T, dif)
    Skpriv=0
    mu=1
    condition = True
    fstruct=jacf
    Tv=lambda x: (np.asmatrix(x)).T
    while condition: #пока не пришли к конвергенции
        Skpriv=Sk
        prevk=k
        Sk=0
        A=np.zeros_like(A)
        b=np.zeros_like(b)

        for i in range(0, len(expdatalist)): #для всех наблюдений
            fstructval=fstruct(expdatalist[i]['x'], k, c, None) #TODO опять же вопрос с неявными функциями
            rt=np.dot (fstructval.T, fstructval)
            A+=rt
            ydif=expdatalist[i]['y']-func(expdatalist[i]['x'],k,c)
            b+=np.dot (fstructval.T, Tv(ydif))   #транспонирование введено для согласования, не коррелирует с формулами
#http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html
        deltak=np.linalg.solve(A,b)  #определяем дельту

        mu=4
        cond2=True
        it=0
        while (cond2):
            Skmu=0
            mu/=2
            for i in range (0, len (expdatalist)):

                vvv=expdatalist[i]['y']-func(expdatalist[i]['x'], mu*deltak.T[0] + k, c)
                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)
                Skmu+=np.dot(vvv.T, vvv)

            it+=1
            if (it>100):
                break
            cond2=Skmu>Skpriv



        k+=mu*deltak.T[0]


        #print (mu, deltak.T[0])

       # k-=0.3*deltak.T[0]

                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)
        Sk=Skmu
        numIterations+=1
        convergence=0

        for i in range (0, M):
            #if (prevk[i]):
            convergence+=math.fabs((deltak[i])/(prevk[i]))

        convergence/=M

        #log+="Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n"


        #print ("Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n")

        if (numIterations>200): #для ради безопасности поставим ограничитель на число итераций
            log+="Break due to max number of iteration exceed"
            print ("GKNU Break due to max number of iteration exceed")
            break
        condition = convergence>math.pow(10, -1*NSIG)
    #print (log)
    #пытаемся проанализировать результат: выводим средний остаток по функции с текущим K
    #по сути тестареа
    # testdiff=0
    #
    # for i in range (0, len(expdatalist)):
    #     testdiff+=math.fabs(func(expdatalist[i]['x'], k,c)[0] - expdatalist[i]['y'][0]) #TODO проверить целесообразность [1]
    # testdiff/=len(expdatalist)
    #
    # print ("testdiff: ", testdiff)

    return k, Sk, numIterations, log




def testNew():
    funstr= ["y[0]+y[1]-y[2]", "y[0]*b[0]-y[1]*b[1]-x[0]-x[1]", "y[1]*b[1]+y[2]*b[2]+x[1]"]
    btrue=[60,40,60]

    updfunstr=list(map(lambda x: x.replace('[','').replace(']',''),  funstr))


    dfdy=lambda x,b,c,y: np.array( [ [1, 1, -1],
                                     [b[0], -b[1], 0],
                                     [0, b[1], b[2]]
    ])

    dfdb=lambda x,b,c,y: np.array ([ [0,    0,    0    ],
                                     [y[0],-y[1], 0    ],
                                     [0,    y[1], y[2] ] ])

    #возвращает функцию
    function=lambda x,b,c,y: [y[0]+y[1]-y[2], y[0]*b[0]-y[1]*b[1]-x[0]-x[1], y[1]*b[1]+y[2]*b[2]+x[1]]

    #возвращает структурную матрицу
    jacf=lambda x,b,c,y: np.dot(dfdb(x,b,c,y), np.linalg.inv(dfdy(x,b,c,y))) #g=dy/db=df/db*np.inv(df/dy)

    #возвращает значение y
    funcf=lambda x,b,c: optimize.root(function, [1, 1, 1], args=(x,b,c),method='lm').x








testNew()



#производство якобиана
    # for i in range (len(updfunstr)):
    #     print (sympy.diff(updfunstr[i], 'b0'), sympy.diff(updfunstr[i], 'b1'), sympy.diff(updfunstr[i], 'b2'))
