__author__ = 'reiner'
#песочница для отработки функций программы, аналога planpos

import random
import math

from scipy import optimize
import numpy as np

import derivations as drv


#уравнения Кирхгофа:

b=(60,60,40) #задаём вектор коэффициентов
w=(1,2) #задаём вектор значений источников питания


#http://pythonworld.ru/tipy-dannyx-v-python/vse-o-funkciyax-i-ix-argumentax.html


def strEvaluator (funstr:list, x:list, b:list=[], c:dict={}):
    """
    Возвращает функцию, которая возвращает list результатов векторной функции, подаваемой во входных параметрах
    Сделано только для неявных функций
    funstr - вектор функций в строковом представлении
    x - вектор входных аргументов
    b - вектор коэффициентов
    c - словарь прочих постоянных
    """
    vardict = {'x':x, 'b':b}
    vardict.update(c)

    def evalfunc(yy):
        #yy - вектор выходных аргументов, которые мы находим методом Ньютона
        vdict = vardict
        vdict['y']=yy



        fff=lambda f: eval (f, None, vardict)

        return list(map(fff, funstr))

    return evalfunc #возвращаем функцию

def debilizm (r:list, let:str):
    """
    Из вектора x[...] формирует словарь {x1: .. ,...} для совместимости с drv.Jacobean
    """
    res=dict()
    for i in range(0, len(r)):
        res[let+str(i)] = r[i]
    return res
#хотели callable - получайте!
def ret_callable_jac (funstr:list, x:list, b:list=[], c:dict={}):
    def innerj (y):
        """

        """
        argseq=list()
        for i in range(0, len(y)):
            argseq.append('y{0}'.format(i))

        updfunstr=list(map(lambda x: x.replace('[','').replace(']',''),  funstr))
        xdict = debilizm(x,'x')
        bdict = debilizm(b,'b')
        cdict = debilizm(c,'c')
        ydict = debilizm(y,'y')
        vardict = xdict
        vardict.update(bdict)
        vardict.update(cdict)
        vardict.update(c)
        vardict.update(ydict)

        return drv.Jakobean (updfunstr, argseq, vardict)

    return innerj

def ret_callable_jac_ultra (funstr:list, x:list, vectorletter:str,  b:list=[], c:dict={}):
    """
    vectorletter - буква вектора, по которому мы берём производную
    """

    def innerj (y):
        """
        y - значение вектора, по которому берём производную

        """
        argseq=list()
        for i in range(0, len(y)):
            argseq.append('{0}{1}'.format(vectorletter, i))

        updfunstr=list(map(lambda x: x.replace('[','').replace(']',''),  funstr))
        xdict = debilizm(x,'x')
        bdict = debilizm(b,'b')
        cdict = debilizm(c,'c')
        ydict = debilizm(y,'y')
        vardict = xdict
        vardict.update(bdict)
        vardict.update(cdict)
        vardict.update(c)
        vardict.update(ydict)

        return drv.Jakobean (updfunstr, argseq, vardict)

    return innerj

def rety (funstr, x, b, c):
    """
    Возвращает y для указанных аргументов
    """
    function = strEvaluator(funstr,x,b,c)
    sol = optimize.root(function, [1, 1, 1], method='lm', jac=ret_callable_jac(funstr, x,b,c))
    return sol.x

def generate_uniform_plan_exp_data(func, xdiapdictlist:list, b:list, c:dict, ydisps:list=None, n=1, outfilename="", listOfOutvars=None):
    """
    Моделирует набор экспериментальных данных, получаемых по равномерному априорному плану.
    Возвращаемое значение - список словарей с ключами x, y, b, c
    Parameters:
    func - callable function, векторная функция, возвращает вектор y при подаче на вход векторов x, b, и словаря с
    xdiapdictlist - вектор диапазонов, заданных словарями 'start 'end'
    b - вектор коэффициентов
    c - словарь дополнительных переменных
    ydisps - вектор дисперсий
    n - объём выборки
    outfilename - имя выходного файла (пока не используется)
    listOfOutvars - список переменных, подлежащих выводу в файл

    example:
        funstr= ["y[0]+y[1]-y[2]", "y[0]*b[0]-y[1]*b[1]-x[0]-x[1]", "y[1]*b[1]+y[2]*b[2]+x[1]"]
    b=[100,200, 300]
    c={}
    for x in  generate_uniform_plan_exp_data(funstr, [{'start':10, 'end':20},{'start':40, 'end':60}], b, c, [0.000001,0, 0.000001], 10):
        print (x)

    """
    for xdiap in xdiapdictlist:  #нашли step
        xdiap['step']=(xdiap['start']-xdiap['end'])/n

        res=list()
    for i in range(0, n):
        appdict=dict()
        xx=list(map (lambda x: x['start']+x['step']*i, xdiapdictlist))

#        print (xx,b,c, funcstrlist)

        yy=func(xx, b)

        #Внесём возмущения:
        if not ydisps==None:
            #distort = lambda i:

            for i in range (0, len(yy)):
                yy[i]=random.normalvariate(yy[i], math.sqrt(ydisps[i]))

            #random.normalvariate(eval(funcstrdict[func], rec), math.sqrt(yvectordispsdict[func])) #иначе просто добавляется одно значение с разбросом

        appdict={'x':xx, 'y':yy}
        res.append(appdict)

    return res

def grandCountGN(funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, vrslst, NSIG=3, kinit=None):
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


    if kinit==None:
        k=np.array((range (1, len(coeffstrlist)+1  ))   )
    else:
        k=kinit


    prevk=k #предыдущее значение вектора коэфф
    convergence=0
    numIterations=1



    A=np.zeros ((len(coeffstrlist), len(coeffstrlist)))
    b=np.zeros((len(coeffstrlist), 1))

    Sk=0
    Skmu=0
    N=len(Xs)  #размер выборки

    ind=0
    func=lambda x,k: np.array(countfunctvect (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, x.tolist(), k.tolist())) #собственно, функция
    for xx in Xs:
        dif=Ys[ind]-np.array(func(xx,k))
        Sk+= np.dot(dif.T, dif)
        ind+=1



    Skpriv=0
    mu=1

    condition = True
    fstruct = lambda x,k: der.Jakobeand (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, x.tolist(), k.tolist())

    Tv=lambda x: (np.asmatrix(x)).T


    while condition: #пока не пришли к конвергенции
        Skpriv=Sk
        prevk=k
        Sk=0
        A=np.zeros_like(A)
        b=np.zeros_like(b)




        for i in range (0, len(Xs)):   #для всех наблюдений
            fstructval=fstruct(Xs[i], k)
            A+=np.dot (fstructval.T, fstructval)
            ydif=Ys[i]-func(Xs[i],k)
            b+=np.dot (fstructval.T, Tv(ydif))   #транспонирование введено для согласования, не коррелирует с формулами

#http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html
        deltak=np.linalg.solve(A,b)  #определяем дельту

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
        testdiff+=math.fabs(func(Xs[i], k)[1] - Ys[i][1
        ])
    testdiff/=len(Xs)


    print ("testdiff: ", testdiff)


    return k, Sk, numIterations, testdiff



   # nvars = len(expdata[0])-len(funcstrstrlist) #количество входных переменных

    #print (nvars)

def retAdvStructJac (funstr:list, x:list, b:list, c:dict={}, yy:list=None):
    """
    Возвращает структурную матрицу для неявной функции
    G=dY/dB=dF/dB * (dF/dY)^-1
    """
    if yy==None:
        y= rety(funstr, x, b, c)
    else:
        y=yy


    #находим dF/dB


    argseq=list()
    for i in range(0, len(b)):
        argseq.append('b{0}'.format(i))

    updfunstr=list(map(lambda x: x.replace('[','').replace(']',''),  funstr))


    xdict = debilizm(x,'x')
    bdict = debilizm(b,'b')
    cdict = debilizm(c,'c')
    ydict = debilizm(y,'y')
    vardict = xdict
    vardict.update(bdict)
    vardict.update(cdict)
    vardict.update(c)
    vardict.update(ydict)



    dFdB=drv.Jakobean (updfunstr, argseq, vardict)

    dFdY=ret_callable_jac (funstr, x, b, c)(y)

    return np.dot(dFdB, np.linalg.inv(dFdY))

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

    # ind=0
    # func=lambda x,k: np.array(countfunctvect (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, x.tolist(), k.tolist())) #собственно, функция
    #
    # for xx in Xs:
    #     dif=Ys[ind]-np.array(func(xx,k))
    #     Sk+= np.dot(dif.T, dif)
    #     ind+=1

    Skpriv=0
    mu=1

    condition = True
#    fstruct = lambda x,k: der.Jakobeand (funcstrdict, invarstrlist, outvarstrlist, coeffstrlist, x.tolist(), k.tolist())
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


def grandCountGN_UltraX (funcf, jacf,  expdatalist:list, kinit:list, NSIG=3):
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
            PP+=np.dot(jacf(expdatalist[i]['x'], k, None), np.transpose(jacf(expdatalist[i]['x'], k, None))  )
            dif = np.array(expdatalist[i]['y'])-np.array(funcf(expdatalist[i]['x'],k))
            Sk+= np.dot(dif.T, dif)
            PYY+=np.dot(jacf(expdatalist[i]['x'], k, None), np.transpose(np.asmatrix(dif)))
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


def test4(): #тестируем на обычной, не неявной функции
    #y1=x0*(b1+b2)
    #y2=x0*b2

    b=[30,40]
    c={}
    #expdata=generate_uniform_plan_exp_data(funstr, [{'start':10, 'end':100},{'start':200, 'end':250}], b, c, [0.000001,0, 0.000001], 50)

    funcf=lambda x,b: np.array ([x[0]*(b[0]+b[1]), x[0]*b[1]])
    jacf = lambda x,b, y: np.matrix([ [x[0], x[0]], [0, x[0] ]   ])

    expdata=generate_uniform_plan_exp_data(funcf, [{'start':10, 'end':100},{'start':200, 'end':250}], b, c, None, 50)

    print  (grandCountGN_Ultra(funcf, jacf, expdata, [1,1]))
    #Почему-то так UltraX c mu выдаёт неверный результат, без mu выдаёт верный, но с большим число итераций

test4()



def test2():
    funstr= ["y[0]+y[1]-y[2]", "y[0]*b[0]-y[1]*b[1]-x[0]-x[1]", "y[1]*b[1]+y[2]*b[2]+x[1]"]
    b=[100,200, 300]
    c={}
    for x in  generate_uniform_plan_exp_data(funstr, [{'start':10, 'end':20},{'start':40, 'end':60}], b, c, [0.000001,0, 0.000001], 10):
        print (x)

#test2()

def test1():
    #на этом простом примере видно, что optimize_root выдаёт верные корни при использовании в связке с strEvaluator
    #однако же при использовании его для модели транзистора резистивной есть расхождения с вариантом функции.
    #возможно, надлежит скормить якобиан однако же




    #funstr= ["x[0]+y[0]+b[0]", "x[1]+y[1]+b[1]" ]
    #x=[1,2]
    #y=[30,30]
    #b=[100,200]
    c={}

    #w преобразовано в х
    funstr= ["y[0]+y[1]-y[2]", "y[0]*b[0]-y[1]*b[1]-x[0]-x[1]", "y[1]*b[1]+y[2]*b[2]+x[1]"]
#    w=[10,60] #задаём
 #   x=w
    b=[60,60,40] #задаём вектор коэффициентов
 #   c=[]
    x=[1,2]

 #   print (ret_callable_jac(funstr, x,b,c)([10,20, 30]))




#    exit(0)


    for i in make_exp_plan_2_generator ((10,20), (60,40),  10):
        x=list(i)
        function = strEvaluator(funstr,x,b,c)

        print (x,b,c)

        #sol = optimize.root(function, [1, 1, 1], method='lm', jac=ret_callable_jac(funstr, x,b,c))

        #TODO можно прогнать тестирование разными методами и выбрать лучший для конкретной задачи
        #sol = optimize.root(function, [1, 1, 1], method='hybr')
        #print (i, sol.x, function(sol.x))


#test1()


    #return strEvaluator(funstr,x,b,c)(y)
def test():
    global b


    print ("points in plan | ", " y vector | ", " left side of equation")

    for i in make_exp_plan_2_generator ((10,20), (60,40),  10):
        w=i
        #sol = optimize.root(kirh, [1, 1, 1 ], jac=jac, method=’hybr’)

        fun = kirhWrapper(b)


        sol = optimize.root(fun, [1, 1, 1 ], method='hybr')


        #можно скормить jacobean. Но он должен быть callable с подачей y на вход, очевидно


        print (i, sol.x, fun(sol.x))



        #print (list(i).extend(sol.x).extend (kirh(sol.x)))


        #написать проверялку, а точно корни-то подходящие?


#print (test1())
#test()

