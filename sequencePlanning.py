__author__ = 'vasilev_is'
import ApriorPlanning as ap, math, numpy as np



def grandCountGN_Ultra (funcf, jacf,  expdatalist:list, kinit:list, NSIG=3):
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
        dif = np.array(expdatalist[i]['y'])-np.array(func(expdatalist[i]['x'],k))
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
            fstructval=fstruct(expdatalist[i]['x'], k, None)
            rt=np.dot (fstructval.T, fstructval)
            A+=rt
            ydif=expdatalist[i]['y']-func(expdatalist[i]['x'],k)
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

                vvv=expdatalist[i]['y']-func(expdatalist[i]['x'], mu*deltak.T[0] + k)
                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)
                Skmu+=np.dot(vvv.T, vvv)

            it+=1
            if (it>100):
                break
            cond2=Skmu>Skpriv

        k+=mu*deltak.T[0]


        print (mu, deltak.T[0])

       # k-=0.3*deltak.T[0]

                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)
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
        testdiff+=math.fabs(func(expdatalist[i]['x'], k)[1] - expdatalist[i]['y'][1]) #TODO проверить целесообразность [1]
    testdiff/=len(expdatalist)

    print ("testdiff: ", testdiff)

    return k, Sk, numIterations, testdiff



def grandSeqPlan (xstart:list, xend:list, N:int, btrue:list, binit:list, c, ydisps, jacf, funcf, NSIG=3):
    """
    Осуществляет последовательное планирование и оценку коэффициентов модели
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: количество точек в плане эксперимента
    :param binit: стартовое значение вектора коэффициентов
    :param btrue: истинное значение вектора коэффициентов
    :param c: словарь дополнительных переменных
    :param ydisps: диагональ ковариационной матрицы y
    :param jacf: функция, возвращающая якобиан модели, входящие x,b,c,y
    :param funcf: функция, возвращающая значение модели, входящие x,b,c
    :param NSIG:
    :return: k, число итераций, лог
    """

    startplan =  ap.makeUniformExpPlan(xstart, xend, N)
    measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c, ydisps)




    b=grandCountGN_Ultra(funcf, jacf, measdata, binit, NSIG)




