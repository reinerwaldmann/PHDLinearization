from Ofiura import Ofiura_general as o_g

__author__ = 'vasilev_is'

import random
import math

import numpy as np


"""

makeUniformExpPlan(xstart:list, xend:list, N:int):
    Делает равномерный неслучайный план эксперимента
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param N: Количество точек в плане
    :return: равномерный план эксперимента

makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    Делает случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор
    :param xstart: начало диапазона x (вектор)
    :param xend: конец диапазона x (вектор)
    :param N: Количество точек в плане

makeMeasAccToPlan(func, expplan:list, b:list, c:dict, ydisps:list=None, n=1, outfilename="", listOfOutvars=None):
     Моделирует измерения в соответствии с планом эксперимента.
    :param func: векторная функция, на вход принимает x, на выходе выдаёт значения y
    :param expplan: план эксперимента (список значений вектора x) (список списков)
    :param b: вектор b - коэффициенты
    :param c: словарь c - дополнительные переменные
    :param ydisps: вектор дисперсий y (диагональ ковариационной матрицы)
    :param n: объём выборки y  - в данной версии не используется
    :param outfilename: имя выходного файла, куда писать план - в данной версии не используется
    :param listOfOutvars: список выносимых переменных - в данной версии не используется
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...

countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    Считает значение определителя Vb, ковариационной матрицы коэффициентов, для данного плана эксперимента
    :param expplan: план эксперимента (список списков)
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента

writePlanToFile (plan, filename='plan.txt'):
    Пишет план в файл. Значения пишутся в формате python (вектор [])
    :param plan: план (список)
    :param filename: имя файла
    :return: ничего

readPlanFromFile (filename='plan.txt'):
    Читает план из файла - запись в виде столбца векторов []
    :param filename: имя файла
    :return: план (список векторов)

def makeUniformExpPlan(xstart:list, xend:list, N:int):
    Создаёт равномерный план эксперимента
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане (внимание! это количество измерений каждого вектора, то есть, реальное кол-во будет N^len(xstart))
    :return: равномерный план эксперимента

makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор

makeMeasAccToPlan(func, expplan:list, b:list, c:dict, Ve=[], n=1, outfilename="", listOfOutvars=None):
    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...


countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    :param expplan: план эксперимента
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента

"""


def for_filter (x,lim):
    for val in x['y']:
        if val>lim:
            return False
    return True

def filterList(inmeasdata, lim=1e55):
    filt = lambda x: for_filter(x,lim)
    return list(filter(filt, inmeasdata))






def writePlanToFile (plan, filename='plan.txt'):
    """
     Пишет план в файл. Значения пишутся в формате python (вектор [])
    :param plan: план (список)
    :param filename: имя файла
    :return: ничего
    """
    if not plan:
        return
    with open (filename, 'wt') as outfile:
        for point in plan:
            outfile.write(point.__str__())
            outfile.write('\n')

def readPlanFromFile (filename='plan.txt'):
    """
    Читает план из файла - запись в виде столбца векторов []
    :param filename: имя файла
    :return: план (список векторов)
    """
    with open (filename, 'rt') as infile:
        lines = infile.readlines()
        return list(map(eval, lines))

def makeUniformExpPlan(xstart:list, xend:list, N:int):
    """
    Создаёт равномерный план эксперимента
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане (внимание! это количество измерений каждого вектора, то есть, реальное кол-во будет N^len(xstart))
    :return: равномерный план эксперимента
    """
    res=list()

    xstartnp=np.array(xstart)
    xendnp=np.array(xend)
    xstep = list((xendnp-xstartnp)/N)

    evalstr="import numpy\n\n"
    lststr=""

    for i in range (len(xstart)):
        evalstr+="\t"*i+"for x{0} in numpy.arange(xstart[{0}], xend[{0}], xstep[{0}]):\n".format(i)
        lststr+="x{0}".format(i)+("," if i+1<len(xstart) else "")
    evalstr+="\t"*(i+1)+"res.append(["+lststr+"])"



    #print (evalstr)
    exec(evalstr,  locals()) #исполняет полученную сроку, собсна, формирует список входных переменных

    # for i in range(N):
    #     res.append(list(xstart+i*xstep))

    # if len(xstart)==1:
    #     for i in range (0, len(res)):
    #         res[i]=[res[i],] #костыль для диода и не нужен



    return res

def makeRandomUniformExpPlan(xstart:list, xend:list, N:int):
    """
    :param xstart: начало диапазона x
    :param xend: конец диапазона x
    :param N: Количество точек в плане
    :return: случайный план эксперимента, в котором точки равномерно распределяются в диапазоне в формате списка словарей 'x': вектор
    """
    res=list()
    for i in range(0, N):
        res.append(o_g.uniformVector(xstart, xend))
    return res

def makeMeasAccToPlan(func, expplan:list, b:list, c:dict, Ve=[], n=1, outfilename="", listOfOutvars=None):
    """

    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """
    res = list()

    for i in range(len(expplan)):
        y=func(expplan[i],b,c)
        if y is None: #если функция вернула чушь, то в measdata её не записывать!
            continue
        #Внесём возмущения:
        if Ve is not None and np.linalg.det(Ve)>10e-15:
            ydisps=np.diag(Ve)
            for k in range(len(y)):
                y[k]=random.normalvariate(y[k], math.sqrt(ydisps[k]))
        curdict = {'x':expplan[i], 'y':y}
        #res[i]["y"]=y
        res.append(curdict)
    return res

def makeMeasOneDot(func, xdot, b:list, c:dict, Ve=[]):
    """

    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """


    y=func(xdot,b,c)
    if y is None: #если функция вернула чушь, то в measdata её не записывать!
        return None
    #Внесём возмущения:
    if Ve is not None and np.linalg.det(Ve)>10e-15:
        ydisps=np.diag(Ve)
        for k in range(len(y)):
            y[k]=random.normalvariate(y[k], math.sqrt(ydisps[k]))
    return y
    #res[i]["y"]=y



def makeMeasOneDot_lognorm(func, xdot, b:list, c:dict, Ve=[]):
    """

    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """


    y=func(xdot,b,c)
    if y is None: #если функция вернула чушь, то в measdata её не записывать!
        return None
    #Внесём возмущения:
    if Ve is not None and np.linalg.det(Ve)>10e-15:
        ydisps=np.diag(Ve)
        for k in range(len(y)):
            #y[k]=random.normalvariate(y[k], math.sqrt(ydisps[k]))
            y[k]=math.exp(random.normalvariate(math.log(y[k]), math.sqrt(ydisps[k])))


    return y



def makeMeasAccToPlan_lognorm(func, expplan:list, b:list, c:dict, Ve=None, n=1, outfilename="", listOfOutvars=None):
    """
    :param func: векторная функция
    :param expplan: план эксперимента (список значений вектора x)
    :param b: вектор b
    :param c: вектор c
    :param Ve: ковариационная матрица (np.array)
    :param n: объём выборки y
    :param outfilename: имя выходного файла, куда писать план
    :param listOfOutvars: список выносимых переменных
    :return: список экспериментальных данных в формате списка словарей 'x':..., 'y':...
    """
    res = list()

    for i in range(len(expplan)):
        y=func(expplan[i],b,c)
        if y is None: #если функция вернула чушь, то в measdata её не записывать!
            continue
        #Внесём возмущения:
        if Ve is not None:
            if np.linalg.det(Ve)>10e-15:
                ydisps=np.diag(Ve)
                for k in range(len(y)):
                    if (y[k]<0):
                        y[k]=-1*math.exp(random.normalvariate(math.log(math.fabs(y[k])), math.sqrt(ydisps[k])))
                    else:
                        y[k]=math.exp(random.normalvariate(math.log(y[k]), math.sqrt(ydisps[k])))

        curdict = {'x':expplan[i], 'y':y}
        #res[i]["y"]=y
        res.append(curdict)
    return res






#FIXME функция makeMeasAccToPlan_lognorm и makeMeasAccToPlan работает в общем  неверно, если говорить об ошибках. И лишь в частном случае диагональной
#ковариационной матрицы ошибок таки правильно, да.
#Как починить? Для каждой точки генерируем вектор случайных величин, коррелированных по матрице с матожиданием 0. Функция делается допиливанием
# generrandvals (M, cov_ksi, nvol=1000). К логарифмированному току прибавляем то, что получилось. Обратно экспоненцируем (как запилено в lognorm)



#генерирует список случайных чисел м=0 д=1
def generlist (n):
    r=list()
    for x in range (0,n):
        r.append(random.gauss(0,1))
    return r



def generrandvals (M, cov_ksi, nvol=1000):
    """
    генерирует кортеж случайных последовательностей, коррелированных по матрице и с матожиданием во вх. параметрах
    """
    M_ksi=np.array(M)
    #проверка корректности:

    if (np.linalg.det (cov_ksi) <0 ):
        print ("Определитель матрицы меньше  0")
        return None


    n=len(M_ksi)  #размерность вектора
    U=list() #генерация случайного вектора U


    for x in range (0,nvol):
        U.append(generlist(n))

    A=np.linalg.cholesky(cov_ksi)

    ksi=list()
    for i in range (0, nvol):
        ksi.append( np.dot(A, np.array(U[i]).T)+M_ksi.T)

    #ksi - список значений случайных векторов, значения векторов коррелированы (значение - вектор)

    #эта часть разбрасывает значения аккуратно по спискам, так, что получается список векторов (вектор - значение)

    res=list()
    for i in range (0, n):
        res.append(list())

    for ni in range (0,nvol):
        for i in range (0, n):
            res[i].append (ksi [ni][i])

    XX=np.array(res)
    COVTEST=np.cov(XX, bias=-1)

    MTEST=list(map (np.mean, res))

    return res, COVTEST, MTEST







def countVbForPlan(expplan:list, b:list,  c:dict, Ve, jac, func=None):
    """
    :param expplan: план эксперимента
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :return: значение определителя для данного плана эксперимента
    """
    G=np.zeros((len(b),len(b))) #матрица G

    for point in expplan:
        jj=np.array(jac(point, b, c, func(point,b,c) if func else None))

          #G+=jj*np.linalg.inv(Ve)*jj.T  #

        #G+=np.dot ( np.dot(jj.T, np.linalg.inv(Ve)), jj)  #поправлено в соответствии с формулой

        G+=np.dot(jj.T, jj) #вставление Ve в середину вроде как особо не требуется


    try:
        return np.linalg.inv(G)
    except BaseException as e:
        print('Fatal error in countVbForPlan: ',e)
        print('b vector=',b)
        print('expplan=',expplan)
        print('G=',G)
        return None


def countVbForMeasdata(b:list,  c:dict, Ve, jac, measdata):
    """
    Для неявных функций актуальным является значение y. В некоторых случаях (напр. при последовательном планировании)
    y уже просчитан, и нет нужды его считать снова, задавая функцию.
    :param expplan: план эксперимента
    :param b: b (вектор коэффициентов)
    :param b: b (вектор коэффициентов)
    :param c: словарь доп. параметров
    :param Ve: ковариационная матрица ошибок экспериментов np.array
    :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
    :param measdata: данные измерений
    :return: значение определителя для данного плана эксперимента
    """
    G=np.zeros((len(b),len(b))) #матрица G

    for point in measdata:
        jj=jac(point['x'], b, c, point['y'])
        #G+=jj*np.linalg.inv(Ve)*jj.T
        #G+=np.dot(jj.T, jj)
        G+=np.dot ( np.dot(jj.T, np.linalg.inv(Ve)), jj)  #поправлено в соответствии с формулой

    try:
        return np.linalg.inv(G)
    except BaseException as e:
        print('Fatal error in countVbForMeasdata: ',e)
        print('b vector=',b)
        print('current point=',point)
        print('G=',G)
        exit(0)


