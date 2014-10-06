__author__ = 'vasilev_is'
import math
import copy

import numpy as np

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



        k+=0.3*deltak.T[0]


        #print (mu, deltak.T[0])

       # k-=0.3*deltak.T[0]

                #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
                # (это матрица с одним столбцом)
        Sk=Skmu
        numIterations+=1
        convergence=0

        for i in range (0, M):
            if (prevk[i]):
                convergence+=math.fabs((deltak[i])/(prevk[i]))

        convergence/=M

        #log+="Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n"


        #print ("Iteration: "+ str(numIterations) + "\n" + "Vect K="+str(k)+"\n"+"Sk="+str(Sk)+"\n\n")

        if (numIterations>200): #для ради безопасности поставим ограничитель на число итераций
            log+="Break due to max number of iteration exceed"
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



def getbSeqPlan (xstart:list, xend:list, N:int, btrue:list, binit:list, c, ydisps, jacf, funcf, NSIG=3):
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

    Для переделывания на реальные измерения btrue=None и функции моделирования переводятся на измерительные
    """

    #делаем ковариационную матрицу y
    Vearr = np.array(ydisps)
    Ve = np.diag(Vearr)


    startplan =  ap.makeUniformExpPlan(xstart, xend, N)
    measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c, ydisps)




    log=""

    b=binit
    for numiter in range(1000000): #ограничитель цикла - если выход произошёл по ограничению, значит, возможна ошибка


        est=grandCountGN_Ultra(funcf, jacf, measdata, b, c, NSIG) #получили оценку b binit=b
        b=est[0]


        Vb=ap.countVbForPlan(startplan, b,  c, Ve, jacf, funcf) #TODO: проработать вопрос с неявными функциями в данном случае!
        #Возможно, надлежит сделать эту функцию универсальной - принимающей и план эксперимента, и measdata.
        #посчитали определитель

        detVb=np.linalg.det(Vb)

        print ("iteration: {0}\nb={1}\ndetVb={2}\nest={3}".format(numiter, b, detVb, (est[1], est[2], est[3])))

        if (detVb<math.pow(10, -6)): #если определитель меньше 10^-6
            return b, numiter, log #то всё вернуть

        #иначе поиск следующей точки плана

        xdot=copy.copy(xstart) #получили начальную точку начальное значение - ровно пополам диапазона
        for i in range(len(xstart)): #присвоили ей значение
            xdot[i]=xstart[i]+(xend[i]-xstart[i])/2

        #объектная функция
        function = lambda x: np.linalg.det(ap.countVbForPlan(ap.appendToList(startplan, x),b,c,Ve,jacf, funcf))
        #каждый раз будет пытаться добавить в план точку и вернуть определитель с добавленной точкой

        xdot=ap.doublesearch (xstart, xend, xdot, function) #оптимизировали значение точки

        startplan.append(xdot)
        measdata.append({'x':xdot, 'y':funcf(xdot,b,c)})


    #окончание этого цикла "естественным путём" говорит о том, что превышено максимальное число итераций
    return b, 1000000, log+"ERROR: maximum number of iteration archieved" #то всё вернуть








def test():
    """
    Тестирует априорное планирование
    :return:
    """
    xstart=[1, 100]
    xend=[20,200]
    N=50
    c={"a":1000}
    funcf=lambda x,b,c: np.array ( [ b[0]+b[1]*x[0]+b[2]*x[1]+b[3]*x[0]*x[1]+b[4]*x[0]*x[0]+b[5]*x[1]*x[1],   b[6]+b[7]*x[0]+b[8]*x[1]+b[9]*x[0]*x[1]+b[10]*x[0]*x[0]+b[11]*x[1]*x[1] ] )
    jacf = lambda x,b,c,y: np.matrix([ [1, x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1], 0, 0, 0, 0, 0, 0],
                                       [0,0,0,0,0,0,1,x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1]] ])
    Ve=np.array([ [0.1, 0],
                  [0, 0.1]]  )
    bstart=[0.8,0.4,1.4,0.2,0.9,0.3,1.4,0.2,2.1,3.1,4.1,5.1]
    blen=  [0.3,0.2,0.2,0.2,0.2,0.3,0.2,0.2]
    bend=  [1.1,0.6,1.6,0.4,1.1,0.6,1.6,0.4,2.5,3.3,4.6,5.6]

    #print (doublesearch ([1, 0.5], [10,10], [9,9], lambda x: x[0]*x[0]+2*x[1]*x[1]+10)) #тестирование поиска

    binit=[x for x in range(len(bstart))]
    print (getbSeqPlan (xstart, xend, N, bstart, bend , c, np.diagonal(Ve), jacf, funcf, NSIG=3))


    # plan=makeUniformExpPlan(xstart, xend, N)
    # func = lambda x,b,c: [x[0]*b[0]+c["a"], x[1]*b[1]+c["a"], x[2]*b[2]+c["a"]]
    # meas = makeMeasAccToPlan(func, plan,  b, c, [0.0001]*3)
    # for x in meas:
    #     print (x)



test()
