__author__ = 'vasilev_is'
import math
import copy, pickle

import numpy as np

import ApriorPlanning as ap

import random

#http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html

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


def getbSeqPlanUltra (xstart:list, xend:list, N:int, btrue:list, binit:list, c, Ve, jacf, funcf, initplan=None, NSIG=3, smallestdetVb=1e-6):
    """
    Осуществляет последовательное планирование и оценку коэффициентов модели
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: количество точек в плане эксперимента
    :param binit: стартовое значение вектора коэффициентов
    :param btrue: истинное значение вектора коэффициентов
    :param c: словарь дополнительных переменных
    :param Ve:  ковариационная матрица y
    :param jacf: функция, возвращающая якобиан модели, входящие x,b,c,y
    :param funcf: функция, возвращающая значение модели, входящие x,b,c
    :param initplan: - начальный план
    :param NSIG:
    :return: k, число итераций, лог
    Для переделывания на реальные измерения btrue=None и функции моделирования переводятся на измерительные
    """
    startplan =  initplan if initplan else ap.makeUniformExpPlan(xstart, xend, N)
    origplan = copy.copy(startplan)
    measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    log=""
    b=binit
    for numiter in range(100): #ограничитель цикла - если выход произошёл по ограничению, значит, возможна ошибка

        est=grandCountGN_Ultra(funcf, jacf, measdata, b, c, NSIG) #получили оценку b binit=b
        b=est[0]
        Sk=est[1]

        Vb=ap.countVbForPlan(startplan, b,  c, Ve, jacf, funcf) #TODO: проработать вопрос с неявными функциями в данном случае!
        #Возможно, надлежит сделать эту функцию универсальной - принимающей и план эксперимента, и measdata.
        #посчитали определитель

        detVb=np.linalg.det(Vb)

        #print ("iteration: {0}\nb={1}\ndetVb={2}\nGKNU output(sum, numiter, log)={3}".format(numiter, b, detVb, (est[1], est[2], est[3])))

        if detVb>0 and detVb<smallestdetVb: #если определитель меньше 10^-6
            return b, numiter, Sk, startplan, origplan, log #то всё вернуть

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
    return b, 1000000, Sk, startplan, origplan, log+"ERROR: maximum number of iterations archieved" #то всё вернуть








def test():
    """
    Тестирует априорное планирование
    :return:
    """
    xstart=[1, 100]
    xend=[20,200]

    N=10
    c={"a":1000}
    funcf=lambda x,b,c: np.array ( [ b[0]+b[1]*x[0]+b[2]*x[1]+b[3]*x[0]*x[1]+b[4]*x[0]*x[0]+b[5]*x[1]*x[1],   b[6]+b[7]*x[0]+b[8]*x[1]+b[9]*x[0]*x[1]+b[10]*x[0]*x[0]+b[11]*x[1]*x[1] ] )
    jacf = lambda x,b,c,y: np.matrix([ [1, x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1], 0, 0, 0, 0, 0, 0],
                                      [0,0,0,0,0,0,1,x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1]] ])

    Ve=np.array([ [0.1, 0],
                  [0, 0.1]]  )
    btrue=[8,4,2,2,9,3,4,2,2,3,4,5]
    bstart=np.array(btrue)-np.array([2]*len(btrue))
    bend=np.array(btrue)+np.array([2]*len(btrue))

    #bstart=[0.8,0.4,1.4,0.2,0.9,0.3,1.4,0.2,2.1,3.1,4.1,5.1]
#    blen=  [0.3,0.2,0.2,0.2,0.2,0.3,0.2,0.2]
    #bend=  [1.1,0.6,1.6,0.4,1.1,0.6,1.6,0.4,2.5,3.3,4.6,5.6]
    binit=[1]*len(btrue)

    #проверяем работу метода оценки
    # startplan =  ap.makeUniformExpPlan(xstart, xend, N)
    # measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    # binit=[1]*len(btrue)
    # print("performing normal research:")
    # startplan =  ap.makeUniformExpPlan(xstart, xend, N)
    # measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c,Ve )
    # gknu=grandCountGN_Ultra (funcf, jacf,  measdata, binit, c, NSIG=3)
    # print ("k, Sk, numIterations, log")
    # print (gknu)
    # print ("AvLog, DLog, sigmaLog")
    # print (ap.logTruthness (measdata, gknu[0], Ve,  funcf, c))
    #
    #
    #
    # print ("\n\nperforming aprior planning:")
    # try: #пытаемся открыть архивированные данные
    #     pkl_file = open('oplan.pkl', 'rb')
    #     oplan = pickle.load(pkl_file)
    #     pkl_file.close()
    # except BaseException:
    #     pkl_file  = open('oplan.pkl', 'wb')
    #     oplan=ap.grandApriornPlanning (xstart, xend, N, bstart, bend, c, Ve, jacf, func=None, Ntries=10)[1]
    #     pickle.dump(oplan, pkl_file, -1)
    #     pkl_file.close()
    # measdata = ap.makeMeasAccToPlan(funcf, oplan, btrue, c,Ve )
    # gknu=grandCountGN_Ultra (funcf, jacf,  measdata, binit, c, NSIG=3)
    # print ("k, Sk, numIterations, log")
    # print (gknu)
    # print ("AvLog, DLog, sigmaLog")
    # print (ap.logTruthness (measdata, gknu[0], Ve,  funcf, c))




    print("\n\nperforming sequence plan:")
    seqplanb=getbSeqPlanUltra (xstart, xend, N, btrue, binit, c, Ve, jacf, funcf)

    print ("b, numiter, Sk, startplan, origplan, log")
    print (seqplanb[0],seqplanb[1],seqplanb[2], seqplanb[5]  )



    print ("AvLog, DLog, sigmaLog")
    measdata = ap.makeMeasAccToPlan(funcf, seqplanb[3], btrue, c,Ve )
    print (ap.logTruthness (measdata, seqplanb[0], Ve,  funcf, c))

    for x in seqplanb[3]:
        print (x)
    print ("\noriginal plan")
    for x in seqplanb[4]:
        print (x)




    # print("\n\nperforming enchanced sequence plan:")
    # #начальный план - априорный
    # seqplanb=getbSeqPlanUltra (xstart, xend, N, btrue, binit, c, Ve, jacf, funcf,initplan=oplan)
    # print ("b, numiter, Sk, startplan, origplan, log")
    # print (seqplanb[0],seqplanb[1],seqplanb[2], seqplanb[5]  )
    # print ("AvLog, DLog, sigmaLog")
    # measdata = ap.makeMeasAccToPlan(funcf, seqplanb[3], btrue, c,Ve )
    # print (ap.logTruthness (measdata, seqplanb[0], Ve,  funcf, c))
    #
    # for x in seqplanb[3]:
    #     print (x)
    # print ("\noriginal plan")
    # for x in seqplanb[4]:
    #     print (x)
    #







test()


def test1():
    xstart=[10, 100]
    xend=[50,300]
    N=50
    btrue=[30,40]
    funcf=lambda x,b, c: np.array ([x[0]*(b[0]+b[1]), x[0]*x[1]*b[1]])
    jacf = lambda x,b, c, y: np.matrix([ [x[0], x[0]], [0, x[0]*x[1] ]   ])
    Ve=np.array([ [0.5, 0],
                  [0, 0.5]]  )
    binit=[1]*len(btrue)
    c={}
    print("performing normal research:")
    startplan =  ap.makeUniformExpPlan(xstart, xend, N)
    measdata = ap.makeMeasAccToPlan(funcf, startplan, btrue, c,Ve )
    gknu=grandCountGN_Ultra (funcf, jacf,  measdata, binit, c, NSIG=3)
    print (gknu)
    print("\n\nperforming sequence plan:")

    seqplanb=getbSeqPlan (xstart, xend, N, btrue, binit, c, Ve, jacf, funcf, NSIG=3, smallestdetVb=1e10)

    print ()
    print (seqplanb[0],seqplanb[1], seqplanb[4]  )

    for x in seqplanb[2]:
        print (x)
    print ("\noriginal plan")
    for x in seqplanb[3]:
        print (x)


#test1()

#Зачем функция оптимизации плана эксперимента добавляет те точки, которые уже в плане эксперимента есть? Наиболее опорные что ли?

