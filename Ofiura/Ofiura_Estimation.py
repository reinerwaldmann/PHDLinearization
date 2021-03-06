__author__ = 'reiner'

import math
import copy

import numpy as np

import Ofiura.Ofiura_Qualitat as o_q
import Ofiura.Ofiura_general as o_g


"""
def grandCountGN_UltraX1 (funcf, jacf,  measdata:list, binit:list, c, NSIG=3):
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param c словарь дополнительных постоянных
    :param NSIG=3 точность (кол-во знаков после запятой)
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения


"""

from scipy import optimize

def make_grad (funcf, jacf,  measdata:list):
    """
    Возвращает градиент объектной функции
    """
    def gradfunk (b):
        gr=None
        for point in measdata:
            dif=np.array(point['y'])-np.array(funcf(point['x'],b))
            jac = jacf(point['x'],b,None,point['y'])

            if gr is None:
                gr=np.dot(jac.T, dif)
            else:
                gr+=np.dot(jac.T, dif)

        return  np.ravel(gr)

    return gradfunk



def make_object_func (funcf, jacf,  measdata:list):
    """
    Возвращает объектную функцию
    """
    def skfunk (b):




        Sk=0
        for point in measdata:



            dif=np.array(point['y'])-np.array(funcf(point['x'],np.array(b)))
            Sk+=np.dot(dif.T,dif)
        return Sk
    return skfunk


def methodnamedecorator(f):
    """
    Хитрый декоратор, который добавляет к словарю, возвращаемому оптимизатором, название метода
    :param f:
    :return:
    """
    def tmp(*args, **kwargs):
        f1= f(*args, **kwargs)
        f1['method']=kwargs['method']
        return f1

    return tmp

@methodnamedecorator
def om (*args, **kwargs):
    """
    Прокладка, вызывающая функцию из оптимизатора, для того, чтоб можно было бы к нему прицепить декторатор
    :param args:
    :param kwargs:
    :return:
    """
    return optimize.minimize(*args, **kwargs)


def optimizer (funcf, jacf,  measdata:list, binit:list, c, NSIG=3, implicit=False, verbose=False, maxiter=1000):

    skfunk = make_object_func(funcf, jacf,measdata)
    grad = make_grad (funcf, jacf,  measdata)


    methods = [om(skfunk, binit, method='Nelder-Mead', jac=False),
               om(skfunk, binit, method='Nelder-Mead', jac=False),

              ]


    #powell = optimize.minimize(skfunk, binit, method='Nelder-Mead', jac=False)
    #neldemead = optimize.minimize(skfunk, binit, method='Nelder-Mead', jac=False)

    ret= min (methods, key=lambda x:  x.fun)
    #ret['method'] =  'sdfsdfsdf'
    return ret









def grandCountGN_UltraX1 (funcf, jacf,  measdata:list, binit:list, c, NSIG=3, implicit=False, verbose=False, maxiter=1000):
    """
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
    :param funcf callable функция, параметры по формату x,b,cTrue
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param c словарь дополнительных постоянных
    :param NSIG=3 точность (кол-во знаков после запятой)
    :param implicit True если функция - неявная, иначе false
    :param verbose Если True, то подробно принтить результаты итераций
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения


    РАБОЧАЯ GEPRUFT!
    """
    #sign - если  1, то b=b+deltab*mu, иначе b=b-deltab*mu. При неявной функции надо ставить sign=0
    sign=0 if implicit else 1

    Sklist=list()
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

            if jac is None:
                return None

            #print (G.shape, jac.T.shape, jac.shape)
            G+=np.dot(jac.T,jac)

            try:
                dif=np.array(point['y'])-np.array(funcf(point['x'],b,c))
            except BaseException as e:
                print('grandCountGN_UltraX1: As funcf returned None, method  stops:', e)
                return None

            # print(dif, jac)
            # print('-----')
            #
            # print()
            if B5 is None:
                B5=np.dot(dif, jac)
            else:
                B5+=np.dot(dif,jac)
            Sk+=np.dot(dif.T,dif)
        #print(np.linalg.inv(G), B5[:,0])
        #костыль для диодной задачи
        if hasattr(B5, 'A1'):
            B5=B5.A1
        try:
            deltab=np.dot(np.linalg.inv(G), B5)
        except BaseException as e:
            print('Error in G:', e)
            print('G=',G)
            return None


        #mu counting
        mu=4
        cond2=True
        it=0
        Skmu=0
        while (cond2):
            Skmu=0
            mu/=2
            for point in measdata:
                try:
                    dif=np.array(point['y'])-np.array(funcf(point['x'],b+deltab*mu,c)) if sign else np.array(point['y'])-np.array(funcf(point['x'],b-deltab*mu,c))
                except:
                    continue

                Skmu+=np.dot(dif.T, dif)
            it+=1
            if (it>100):
                log+="Mu counting: break due to max number of iteration exceed"
                break
            cond2=Skmu>Sk

        b=b+deltab*mu if sign else b-deltab*mu

        #Sk=Skmu

        Sklist.append(Sk)
        if verbose:
            print ("Sk:",Sk)
            print ("Iteration {0} mu={1} delta={2} deltamu={3} resb={4}".format(numiter, mu, deltab, deltab*mu, b))

        numiter+=1

        condition=False
        for i in range (len(b)):
            if math.fabs ((b[i]-bpriv[i])/bpriv[i])>math.pow(10,-1*NSIG):
                condition=True


        if numiter>maxiter: #max number of iterations
            log+="GKNUX1: Break due to max number of iteration exceed"
            break


    return b, numiter, log, Sklist, Sk



def selectBestEstim (gknuxdictsarr:list):
    """
    1. Фильтрует массив оценок, выкидывая плохие
    2. Выбирает наилучший, возвращает наилучший (или просто сортированный)
    :param gknuxdictsarr
    """
    #1. Фильтр входного массива по чему-то там
    def sortingfunc (arg):
        return arg['AvLogTruth']

    #2 опции


    resl=sorted(gknuxdictsarr, key=sortingfunc)

    if resl:
        return resl[0]
    else:
        print ('Не удалось адекватно оценить параметры. Попробуйте больше интераций')
        return None
    #eturn gknuxdictsarr.sort(key=sortingfunc)[0]

def grandCountGN_UltraX_ExtraStart (funcf, jacf,  measdata:list, bstart:list, bend:list, c, Ve,  NSIG=3, implicit=False, verbose=False, Ntries=10, name=''):
    """
    Новая реализация оценки, в соответствии с методом выбора начального приближения (22.12.2014 tasklist)
    """
    resarray=list()
    for i in range (Ntries):
        #получение uniform binit
        binit = o_g.uniformVector(bstart, bend)

        try:
            resarray.append(grandCountGN_UltraX_Qualitat (funcf, jacf,  measdata, binit, c, Ve,  NSIG, implicit, verbose, name))
        except BaseException:
            pass


    return selectBestEstim (resarray)


def grandCountGN_UltraX (funcf, jacf,  expdatalist:list, kinit:list, c=None, NSIG=3):
    """
    Подгоняет коэфф. методом Ньютона-Гаусса с изменяемой длиной шага (через mu)
    тестировалась на явной функции когда-то (?)
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


            PP+=np.dot(jacf(expdatalist[i]['x'], k, None).T, (jacf(expdatalist[i]['x'], k, None))  )

            dif = np.array(expdatalist[i]['y'])-np.array(funcf(expdatalist[i]['x'],k))
            Sk+= np.dot(dif.T, dif)
            PYY+=np.dot(jacf(expdatalist[i]['x'], k, None), dif.T)
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
        testdiff+=math.fabs(funcf(expdatalist[i]['x'], k)[1] - expdatalist[i]['y'][1]) #TODO проверить целесообразность [1]
    testdiff/=len(expdatalist)


    print ("testdiff: ", testdiff)


    #return k, Sk, numIterations, testdiff
    return k, numIterations, log




def grandCountGN(funcf, jacf,  measdata:list, binit:list, c, NSIG=3):
    """
   Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param c словарь дополнительных постоянных
    :param NSIG=3 точность (кол-во знаков после запятой)
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
    """
    log=""#строка, куда пишутся всякие сообщения
    k=binit
    prevk=k #предыдущее значение вектора коэфф
    convergence=0
    numIterations=1
    A=np.zeros ((len(k), len(k)))
    b=np.zeros((len(k), 1))
    Sk=0
    Skmu=0
    N=len(measdata)  #размер выборки
    ind=0
    for xx in measdata:
        dif=measdata[ind]['y']-np.array(funcf(xx,k,c))
        Sk+= np.dot(dif.T, dif)
        ind+=1
    Skpriv=0
    mu=2
    condition = True

    Tv=lambda x: (np.asmatrix(x)).T

    while condition: #пока не пришли к конвергенции
        Skpriv=Sk
        prevk=k
        Sk=0
        A=np.zeros_like(A)
        b=np.zeros_like(b)

        for i in range (len(measdata)):   #для всех наблюдений
            fstructval=jacf(measdata[i]['x'], k, measdata[i]['y']) #значение якобиана в данной точке

            A+=np.dot (fstructval.T, fstructval)

            dif=np.array(measdata[ind]['y'])-np.array(funcf(xx,k,c))
            #b+=np.dot (fstructval.T, Tv(dif))   #транспонирование введено для согласования, не коррелирует с формулами
            b+=np.dot(dif,fstructval) #по формуле из gknux1  #dif - это вектор, возможно, в итоге придётся нечто транспонировать

#http://docs.scipy.org/doc/numpy/reference/generated/numpy.linalg.solve.html
        deltak=np.linalg.solve(A,b)  #определяем дельту

        mu=0.5

        # cond2=True
        # it=0
        # while (cond2):
        #     Skmu=0
        #     mu/=2
        #     for i in range (0, len (Xs)):
        #
        #         vvv=Ys[i]-func(Xs[i], mu*deltak.T[0] + k)
        #         #почему так? потому, что numpy.linalg.solve выдаёт вертикальный массив, трактуемый как список списков
        #         # (это матрица с одним столбцом)
        #
        #
        #         Skmu+=np.dot(vvv.T, vvv)
        #
        #
        #     it+=1
        #     if (it>1000):
        #         break
        #     cond2=Skmu>Skpriv

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

