__author__ = 'reiner'

import math
import copy

import mpmath as mpm

#Finished


#во всём проекте для матриц mpmath примем такую точность
mpm.dps = 40
mpm.pretty = True



def makeSkInitMpmath (funcf, measdata, binit, c):
    Sk=mpm.mpf(0.0)#или mpm.mpf(0.0)
    for point in measdata:
        fxbc=funcf(point['x'],binit,c)
        if fxbc is None:
            print ("ERROR: makeSkInit: funcf returned None")
            return None
        dif=point['y']-fxbc
        Sk+=(dif.T*dif)[0]
    return Sk


def makeAinitMpmath (bstart, bend, Skbinit, binit, isBinitGood=True):
    """
    расчёт ведётся по правилу "binit есть хорошее предположение в рамках области"
    :param bstart:
    :param bend:
    :param Skbinit:
    :param binit:
    :return:
    """
    A=mpm.matrix( len(binit), 2 )
    for i in range (len(binit)):
        if isBinitGood:
            A[i][0] = mpm.mpf('0.001')*(binit[i]-bstart[i])*Skbinit
            A[i][1] = mpm.mpf('0.001')*(bend[i]-binit[i])*Skbinit
        else:
            A[i][0] = A[i][1] = mpm.mpf('0.001')*(bend[i]-bstart[i])*Skbinit  #универсально
    return A

def countSklimsMpmath(A,b,bstart, bend):
    partone=parttwo=mpm.mpf(0)
    for i in range (len(b)):
        partone+=A[i][0]/(b[i]-bstart[i]) #UNSAFE: если вдруг, ну вдруг b=bstart или bend, то будет <strike> треш </strike> креш с делением на ноль.
        parttwo+=A[i][1]/(bend[i]-b[i])
    return partone+parttwo

def countNMpmath (A, b, bstart, bend):
    N=mpm.matrix((len(b),len(b)))
    for j in range (len(b)):
        partone=2*A[j][0]/(b[j]-bstart[j])**3
        parttwo=2*A[j][1]/(bend[j]-b[j])**3
        N[j][j]+=parttwo+partone #так как матрица нулевая
    return N



def  grandCountGN_UltraX1_Limited_wrapperMpmath (funcf, jacf,  measdata:list, binit:list, bstart:list, bend:list, c, NSIG=50, NSIGGENERAL=50, implicit=False, verbose=False, verbose_wrapper=False, isBinitGood=True):
    """
    Обёртка для grandCountGN_UltraX1_Limited для реализации общего алгоритма
    :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param bstart:list нижняя граница b
    :param bend:list верхняя граница b
    :param c словарь дополнительных постоянных
    :param A матрица коэффициентов a
    :param NSIG=3 точность (кол-во знаков после запятой)
    :param implicit True если функция - неявная, иначе false
    :param verbose Если True, то подробно принтить результаты итераций
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
    """

    maxiter=100
    b,bpriv=binit,binit
    gknux=None
    gknuxlist=list()
    Skinit = makeSkInitMpmath(funcf,measdata,binit, c)
    A=makeAinitMpmath(bstart, bend,Skinit,binit, isBinitGood)
    log=''

    if verbose_wrapper:
        print ('==grandCountGN_UltraX1_Limited_wrapperMpmath is launched==\n\n')

    for numiter in range (maxiter):
        bpriv=copy.copy(b)
        gknux=grandCountGN_UltraX1_Limited_wrapperMpmath (funcf, jacf,  measdata, b, bstart, bend, c, A, NSIG, implicit, verbose) #посчитали b
        if gknux is None:
            print ("grandCountGN_UltraX1_Limited_wrapper crashed on some iteration")
            continue
        gknuxlist.append(gknux)
        if verbose_wrapper:
            print ('Iteration \n',numiter,'\n' ,gknux)
        b=gknux[0]
        if not gknux[2]=='':
            #log+="On gknux iteration "+numiter+": "+ gknux[2]+"\n"
            log+="On gknux iteration {0}: {1}\n".format (numiter, gknux[2])
        for j in range (len(binit)): #уменьшили в два раза
            A[j][0]*=mpm.mpf('0.5')
            A[j][1]*=mpm.mpf('0.5')
        condition=False
        for i in range (len(b)):
            if mpm.fabs ((b[i]-bpriv[i])/bpriv[i]) > math.pow(10,-1*NSIGGENERAL):
                condition=True
        if not condition:
            break
        #мол если хоть один компонент вектора b значимо изменился, тогда продолжать. Иначе программа дойдёт до break и цикл прекратится

#    print ('grandCountGN_UltraX1_Limited_wrapper iterations number:', numiter)
#    print (gknux[0], gknux[1], log, gknux[3], gknux[4])
    return gknux[0], gknux[1], log, gknux[3], gknux[4]



def grandCountGN_UltraX1_mpmath_Limited (funcf, jacf,  measdata:list, binit:list, bstart:list, bend:list, c, A, NSIG=3, implicit=False, verbose=False):
    """
    Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
    В стандартный поток вывода выводит отладочную информацию по каждой итерации
     :param funcf callable функция, параметры по формату x,b,c
    :param jacf callable функция, параметры по формату x,b,c,y
    :param measdata:list список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
    :param binit:list начальное приближение b
    :param bstart:list нижняя граница b
    :param bend:list верхняя граница b
    :param c словарь дополнительных постоянных
    :param A матрица коэффициентов a
    :param NSIG=3 точность (кол-во знаков после запятой)
    :param implicit True если функция - неявная, иначе false
    :param verbose Если True, то подробно принтить результаты итераций
    :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения

    """

    Sklist=list()
    b=binit
    log=""
    numiter=0
    condition=True

    def throwError (msg):
        #global b, numiter, log, Sklist, Sk
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
            G+=jac.T*jac

            fxbc=funcf(point['x'],b,c)
            if fxbc is None:
                throwError ("Funcf is None")
            dif=point['y']-fxbc

            if B5 is None:
            #if B5==mpm.matrix(1,m):
                B5=dif*jac
            else:
                B5+=dif*jac
            Sk+=(dif.T*dif)[0]

            #Вставка для работы с ограничениями
            N = countNMpmath(A,b, bstart, bend)
            Sklims = countSklimsMpmath(A,b,bstart,bend)
            G=G+N
            Sk+=Sklims

        try:
            deltab=(G**-1)*B5.T
        except BaseException as e:
            print('G=',G)
            print('B5=',B5)
            throwError('Error in G:'+e.__str__())

        #mu counting
        mu=mpm.mpf(4)
        cond2=True
        it=0
        while (cond2):
            #Skmu=mpm.mpf(0)
            Skmu=countSklimsMpmath(A,b,bstart,bend)
            mu/=mpm.mpf(2)
            for point in measdata:
                try:
                    dif=point['y']-funcf(point['x'],b-deltab*mu,c) if implicit else point['y']-funcf(point['x'],b+deltab*mu,c)
                except:
                    log+='\nmu counting dif counting throwed an exception while {0}'.format((point, b, b-deltab*mu))
                    continue
                    #весьма спорный вариант, именно из-за него, мне кажется, и бывает mu counting iteration exceed
                    #типа если точка не считается, то отклонить эту точку
                Skmu+=(dif.T*dif)[0]
            it+=1
            if (it>100):
                log+="Mu counting: break due to max number of iteration exceed"
                break
            cond2=Skmu>Sk
        b=b-deltab*mu if implicit else b+deltab*mu
        Sklist.append(Sk)
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
