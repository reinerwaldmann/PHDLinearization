__author__ = 'SW'
import math
import copy

import mpmath as mpm


#во всём проекте для матриц mpmath примем такую точность
mpm.dps = 50;
mpm.pretty = True

def grandCountGN_UltraX1_mpmath (funcf, jacf,  measdata:list, binit:list, c, NSIG=3, implicit=False, verbose=False):
    """
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
    :param funcf callable функция, параметры по формату x,b,c, на выходе вектор (pure Python)
    :param jacf callable функция, параметры по формату x,b,c,y, на выходе матрица mpmath.matrix
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param c словарь дополнительных постоянных
    :param NSIG=3 точность (кол-во знаков после запятой)
    :param sign - если  1, то b=b+deltab*mu, иначе b=b-deltab*mu. При неявной функции надо ставить sign=0
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
    РАБОЧАЯ GEPRUFT!
    """

    Sklist=list()
    b=binit
    log=""
    numiter=0
    condition=True

    def throwError (msg):
        global b, numiter, log, Sklist, Sk
        log+=msg
        return b, numiter, log, Sklist, Sk

    while (condition):
        m=len(b) #число коэффициентов
        #G=np.zeros((m,m))
        G=mpm.matrix(m)
        #B5=mpm.matrix(1,m)
        B5=None
        bpriv=copy.copy(b)
        Sk=mpm.mpf(0)
        for point in measdata:
            jac=jacf(point['x'],b,c,point['y'])

            if jac is None:
                 throwError ("Jac is None")

            #G+=np.dot(jac.T,jac)
            G+=jac.T*jac
            #dif=np.array(point['y'])-np.array(funcf(point['x'],b,c))
            fxbc=funcf(point['x'],b,c)

            if fxbc is None:
                throwError ("Jac is None")
            dif=point['y']-fxbc

            if B5 is None:
            #if B5==mpm.matrix(1,m):
                B5=dif*jac
            else:
                B5+=dif*jac
            Sk+=dif.T*dif
        try:
            deltab=(G**-1)*B5
        except BaseException as e:
            print('G=',G)
            print('B5=',B5)
            throwError('Error in G:', e)

        #mu counting
        mu=mpm.mpf(4)
        cond2=True
        it=0
        while (cond2):
            Skmu=mpm.mpf(0)
            mu/=mpm.mpf(2)
            for point in measdata:
                dif=point['y']-funcf(point['x'],b-deltab*mu,c) if implicit else point['y']-funcf(point['x'],b+deltab*mu,c)
                Skmu+=dif.T*dif
            it+=1
            if (it>100):
                log+="Mu counting: break due to max number of iteration exceed"
                break
            cond2=Skmu>Sk
        b=b-deltab*mu if implicit else b+deltab*mu

        if verbose:
            print ("Sk:",Sk)
            print ("Iteration {0} mu={1} delta={2} deltamu={3} resb={4}".format(numiter, mu, deltab, deltab*mu, b))

        numiter+=1

        condition=False
        for i in range (len(b)):
            if mpm.fabs ((b[i]-bpriv[i])/bpriv[i])>math.pow(10,-1*NSIG):
                condition=True

        if numiter>500: #max number of iterations
            throwError("Break due to max number of iteration exceed")

    return b, numiter, log, Sklist, Sk
