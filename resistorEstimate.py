

__author__ = 'vasilev_is'

"""
Назначение файла - тестирование возможности оценки параметров транзистора с помощью метода Ньютона-Гаусса
Функции - метод Ньютона-Гаусса в исходном виде, функция, якобиан, априорное планирование

Сперва возьмём обычную резистивную модель и попробуем её оценить.  С помощью обычного (равномерного) планирования, с помощью априорного плана
"""

import math
import copy

import numpy as np
from scipy import optimize

import ApriorPlanning as ap


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
            fstructval=fstruct(expdatalist[i]['x'], k, c, expdatalist[i]['y'])
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




        print (mu, deltak.T[0], k , Sk)

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

        if (numIterations>100): #для ради безопасности поставим ограничитель на число итераций
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



def grandCountGN_UltraX1 (funcf, jacf,  measdata:list, binit:list, c, NSIG=3):
    """
    WORKING: PROVED!

    """
    #TODO: добавить mu, ибо расхождение при Ve!=None

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
            G+=np.dot(jac.T,jac)
            dif=point['y']-funcf(point['x'],b,c)
            if B5==None:
                #B5=np.dot(jac,dif)
                B5=np.dot(dif, jac)
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




def grandCountGN_UltraX (funcf, jacf,  expdatalist:list, kinit:list, c, NSIG=3):
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
            PP+=np.dot(jacf(expdatalist[i]['x'], k, c, expdatalist[i]['y']), np.transpose(jacf(expdatalist[i]['x'], k, c, expdatalist[i]['y']))  )
            dif = np.array(expdatalist[i]['y'])-np.array(funcf(expdatalist[i]['x'],k,c))
            Sk+= np.dot(dif.T, dif)
            PYY+=np.dot(jacf(expdatalist[i]['x'], k, c, expdatalist[i]['y']), np.transpose(np.asmatrix(dif)))

        print (PP)

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


        #for i in range (0, len(k)):
        #    k[i]+=deltak.transpose()[0][i]*mu


        print ('mu', mu)


        Sk=Skmu

    #***************





        #k+=mu*deltak.T[0]
        print (k)



        #k-=0.3*deltak.T[0]


                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)




        #Sk=Skmu


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


    return k, Sk, numIterations, testdiff



def jjacf (x,b,c,y, dfdb, dfdy):

    q=dfdy(y,x,b,c)
    k=dfdy(y,x,b,c).shape[0]
    # for i in range(k):
    #     for j in range(k):
    #         d=q[i,j]
    #         q[j,i]=d
    #сделано как у Попова, но теоретически подобное траспонирование ошибочно!!!!

    q=np.linalg.inv(q)
    return np.dot(q,  dfdb(y,x,b,c))


def testNew():
    funstr= ["y[0]+y[1]-y[2]", "y[0]*b[0]-y[1]*b[1]-x[0]-x[1]", "y[1]*b[1]+y[2]*b[2]+x[1]"]


    updfunstr=list(map(lambda x: x.replace('[','').replace(']',''),  funstr))


    dfdy=lambda y,x,b,c: np.array( [ [1, 1, -1],
                                     [b[0], -b[1], 0],
                                     [0, b[1], b[2]]
    ])

    dfdb=lambda y,x,b,c: np.array ([[ 0,    0,    0    ],
                                     [y[0],-y[1], 0    ],
                                     [0,    y[1], y[2] ] ])

    #возвращает функцию
    function=lambda y,x,b,c: [y[0]+y[1]-y[2], y[0]*b[0]-y[1]*b[1]-x[0]-x[1], y[1]*b[1]+y[2]*b[2]+x[1]]

    #возвращает структурную матрицу
    jacf=lambda x,b,c,y: jjacf(x,b,c,y,dfdb,dfdy)

    #jacf=lambda x,b,c,y: np.dot(dfdb(y,x,b,c), np.linalg.inv(dfdy(y,x,b,c))) #g=dy/db=df/db*np.inv(df/dy)



    #возвращает значение y
    funcf=lambda x,b,c: optimize.root(function, [1, 1, 1], args=(x,b,c),method='lm', jac=dfdy).x


    #теперь попробуем сделать эксперимент.
    c={}
    Ve=np.array([ [0.00001, 0, 0],
                  [0, 0.00001, 0],
                  [0, 0, 0.00001]  ]  )

    btrue=[60,60,40]
    bstart=np.array(btrue)-np.array([2]*len(btrue))
    bend=np.array(btrue)+np.array([2]*len(btrue))
    binit=[60,55,45]

    xstart=[10,40]
    xend=[20,60]

    N=10

    print("performing normal research:")
    startplan =  ap.makeUniformExpPlan(xstart, xend, N)
    measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)


    gknu=grandCountGN_UltraX1 (funcf, jacf,  measdata, binit, c, NSIG=6)
    print (gknu)
    print (ap.getQualitat(measdata, gknu[0], Ve,  funcf, c))


    N=20
    print ("\n\nperforming aprior plan")
    oplan=ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, func=funcf, Ntries=6)[1]
    ap.writePlanToFile(oplan)
    measdata = ap.makeMeasAccToPlan(funcf, oplan, btrue, c,Ve )
    gknu=grandCountGN_Ultra (funcf, jacf,  measdata, binit, c, NSIG=3)
    print (gknu)
    print (ap.getQualitat(measdata, gknu[0], Ve,  funcf, c))


testNew()



#производство якобиана
    # for i in range (len(updfunstr)):
    #     print (sympy.diff(updfunstr[i], 'b0'), sympy.diff(updfunstr[i], 'b1'), sympy.diff(updfunstr[i], 'b2'))
