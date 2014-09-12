from copy import deepcopy

__author__ = 'vasilev_is'


import numpy as np
import random
import math
"""
В этом файле реализованы функции генерации случайных значений векторных функций.
"""


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



def generate (funcstrdict, xvectorlistsdict, spreadvarslist=None, Vx=None, nvolx=1, yvectordispsdict=None, nvoly=1, outfilename="", listoutvars=None):
    """
Функция способна в широких пределах моделировать экспериментальные данные, получаемые в процессе исследования
сущности, задаваемой многооткликовой моделью.
Описание прилагается

Входные параметры:
словарь входных функций  - выходная переменная к строке функции
словарь значений входных переменных  - входная переменная к списку значений (возможно, задаваемому как arange, range или каким-либо ещё образом). Если вх. перем. имеют разброс, то это список матожиданий
список переменных с разбросом, в строковом представлении
V, ковариационная матрица переменных с разбросом. Диагональная V говорит о наличии дисперсии, но отсутствии ковариации
nvolx, размер выборки входных переменных, получаемых с разбросом, для каждого матожидания
словарь дисперсий выходных переменных - выходная переменная - дисперсия.
nvoly, размер выборки выходных переменых, получаемых с разбросом, для каждого матожидания
имя выходного файла, если имеется потребность в записи входных и выходных переменных в выходной файл
listoutvars  -  список переменных, подлежащих выводу в файл. Если None, выводится всё


Выходные значения:
список словарей входная/выходная переменная - значение
опционально файл, где записаны входные переменные - выходные переменные. Столбцы озаглавливаются.



#TODO:  Есть режим, при котором необходимо эмулировать ошибку (погрешность) аппаратуры напрямую.
Ошибка в данном случае аддитивна. (то есть Урез = Уисх+Е, где Е определяется исходя из погрешности прибора.
Вопрос: 1. Верно ли, что получаемые этой функцией в этом виде данные измерений (то есть, выборка нормального распределения
с дисперсией задаваемой и средним истинным значением измерения) эмулируют ошибку аппаратуры верно
2. Если так, то как показатели точности аппаратуры переводятся в значений оной дисперсии?
3. Если задана относительная ошибка измерений, скажем, точность сост. 10%, то как это описать дисперсией (и возможно ли?)
Скорее всего, требуется иной подход, именно аддитивный - прибавление к точной части некоторой ошибки, величиная которой
случайна, знак - тоже, но она не превышает 10% от исходных точных данных. 
    """



    evalstr=""
    counter=0
    for key in xvectorlistsdict.keys():
        for i in range (0, counter):
            evalstr+="\t"
        evalstr+="for "+key+" in "+xvectorlistsdict[key].__str__()+":\n"
        counter+=1
    #построена строка с циклами

    varslist=list()

    for i in range (0, counter):
        evalstr+="\t"
    evalstr+="resdict=dict()\n"
    for key in xvectorlistsdict.keys():
        for i in range (0, counter):
            evalstr+="\t"
        evalstr+="resdict[\""+key+"\"]="+key
        evalstr+="\n"
    for i in range (0, counter):
        evalstr+="\t"
    evalstr+="varslist.append(resdict)"
#    print (evalstr)

    #    resdict[key]=list()
     #   evalstr+="resdict[\""+key+"\"].append("+funcstrdict[key]+")\n"

    exec(evalstr, xvectorlistsdict, locals()) #исполняет полученную сроку, собсна, формирует список входных переменных

    #print (varslist)
    #varslist  - список словарей входных переменных

     #сюда должен быть добавлен вариант с генерацией рандомных последователььностей входных переменных
     #есть функция, которая генерирует рандомные коррелированные последовательности.

    #если задан разброс входных параметров
    #то мы вот берём и на каждое значение входных параметров генерим их выборку со значением в качестве матожидания, Vx в качестве корреляционной матрицы

    rndinvalslist=list() #сначала всё новые значения сложить в этот список, потом  уже добавлять в varslist, иначе бесконечность
    if spreadvarslist!=None and  Vx!=None:
        for rec in varslist: #для каждой записи из списка значений входных переменных
            listForM=list()  #составить вектор M
            for var in spreadvarslist: #для каждой переменной из распределённых добавить её значение в вектор М
                listForM.append(rec[var])
            M=np.array(listForM) #сделали аррэй из списка вообще можно б и сразу в аррэй
            rndvals = generrandvals (M, Vx, nvolx)[0] #получили кортеж случайных последовательностей


            for m in range (0, nvolx):  #проходим по каждому случайному вектору
                newrec = deepcopy(rec) #копируем запись
                for ii in range (0, len (spreadvarslist)):#для каждой переменной из распределённых
                    newrec[spreadvarslist[ii]]=rndvals[ii][m]
                    pass


                rndinvalslist.append(newrec)
            varslist.remove(rec) # удаляем  исходное, без разброса

        varslist+=rndinvalslist #добавляем к списку то, что мы получили с разбросом





    #посчитаем значения функции

    for rec in varslist: #для каждой записи из списка значений входных переменных, полученного на предыдущем этапе
        for func in funcstrdict.keys():  #для каждой функции из словаря функций

            if (yvectordispsdict==None):   #если не задано дисперсий для выходных функций (! если уж заданы дисперсии, то должны быть заданы для всех функций!)
                #если надо просто получить y, без разброса
                rec[func]=eval(funcstrdict[func], rec)
            else: # если дисперсии y таки заданы
                if (nvoly>1): #если сказано получать несколько значений y для данных входных перем
                    rec[func]=list()  # то значение y становится списком
                    for i in range(0, nvoly):
                        rec[func].append (random.normalvariate(eval(funcstrdict[func], rec), math.sqrt(yvectordispsdict[func]))) #и в него добавляются значения

                else:
                    rec[func]=random.normalvariate(eval(funcstrdict[func], rec), math.sqrt(yvectordispsdict[func])) #иначе просто добавляется одно значение с разбросом

            del rec["__builtins__"]  #с какого эта запись вообще пишется в словарь???? ОЛОЛО






    # запись в файл


    if len (varslist)==0:
        return varslist


    try:
        with open(outfilename, "wt") as file:

            for key in varslist[0].keys():
                if listoutvars!=None:
                    if key in listoutvars:
                        file.write(key)
                        file.write("\t")
                else:
                    file.write(key)
                    file.write("\t")
            file.write("\n")


            for rec in varslist:
                for key in rec.keys():
                    if listoutvars!=None:
                        if key in listoutvars:
                            file.write(rec[key].__str__())
                            file.write("\t")
                    else:
                        file.write(rec[key].__str__())
                        file.write("\t")
                file.write("\n")
    except BaseException:
        pass



    return varslist


def generateAsDict (funcstrdict, xvectorlistsdict, spreadvarslist=None, Vx=None, nvolx=1, yvectordispsdict=None, nvoly=1, outfilename=""):
    """
        режем список словарей на словарь списков

    """


    vrslst=generate (funcstrdict, xvectorlistsdict, spreadvarslist, Vx, nvolx, yvectordispsdict, nvoly, outfilename)

    if len(vrslst)==0:
        return None

    #режем список словарей на списки, соответствующие идентификаторам словарям - из этого составить функцию и добавить в генератор
    resdict=dict()
    for key in vrslst[0].keys():
        resdict[key]=list()

    for item in vrslst:
        for key in item.keys():
            resdict[key].append(item[key])

    return resdict





#test area:
#
#
# funcstrdict= {"y1":"u1* (r2+r3)/(r1+r2+r3)", "y2":"u1* r3/(r1+r2+r3)"}
# #xvectorlistsdict = {"u1":[100],  "r1":list(np.arange(0,5,1)), "r2":list(range(0,3,1)), "r3":list(range (20,30,5))}
#
# xvectorlistsdict = {"u1":[100],  "r1":[20], "r2":[30], "r3":[400]}
# spreadvarslist  = ["r1", "r2", "r3"]
# V=np.array       ( [[4, 2, 3],
#                     [2, 9, 6],
#                     [3, 6, 16]])
#
#
# #{"y1":16, "y2":25 }
# vrslst=generate (funcstrdict, xvectorlistsdict, spreadvarslist, Vx=V, nvolx=100, yvectordispsdict=None, nvoly=1, outfilename="t.txt")
#
# for item in vrslst:
#     print (item)
#
# M=np.array([20,30,400])
#
#
# #print (generrandvals (M, V, nvol=10)[0] )











#
# def test (inp):
#     for j in list (inp):
#         print (j)


#test (np.arange(0,1,0.1))



