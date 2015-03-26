__author__ = 'vasilev_is'

import pickle

import numpy as np

import Ofiura.Ofiura_general as o_g


#конкретно к модели диода относящиеся вещи

btrue = [1.238e-14, 1.8, 100]
leny=1 #длина выходного вектора

#print (np.linalg.norm(btrue))


def makeplan_lambda (plantype, n, bstart, bend):
    global btrue
    #- делает нужный план в зависимости от параметров plantype и n
    #Надо создать пачку априорных планов - 5 10 15 20 25 30 35 40 значений

    if plantype: #если априорный
        pass
    else:
        pass

def makebinit_lambda (bstart, bend):
    global btrue
    #- делает начальное значение в зависимости от параметров lbinitbtrue
    binit = o_g.uniformVector(bstart, bend)
    return binit, np.linalg.norm(np.array(binit)-np.array(btrue))



def makediap_lambda_absolute (diapwidth, assym):
    #- делает диапазон в зависимости от параметров diapwidth, assym
    # diapwisth, assym задаются абсолютными диапазонами.
    global btrue
    btrue1=np.array(btrue)
    m = len(btrue) #длина вектора b
    centervector=np.array([assym/(m**.5)]*m) + btrue1
    halfrange = np.array([.5*diapwidth/(m**.5)]*m)
    bstart, bend = centervector-halfrange, centervector+halfrange
    return bstart, bend

def makediap_lambda (diapwidth, assym):
    """ делает диапазон в зависимости от параметров diapwidth, assym
     diapwisth, assym задаются в долях от btrue """
        #тестовый сет для этой функции
    # btrue1=np.array(btrue)
    # bassym = .1  #дельта, расстояние, на которое надо отойти от оригинального вектора
    # diapwidth=.40
    # crisp=makediap_lambda(diapwidth, bassym)
    #
    # print (crisp)
    # print ((crisp[1]-crisp[0])/2+crisp[0]-bassym*btrue1)

    global btrue
    btrue1=np.array(btrue)
    m = len(btrue) #длина вектора b
    centervector=btrue1*assym + btrue1
    halfrange = btrue1*diapwidth/2
    bstart, bend = centervector-halfrange, centervector+halfrange
    return bstart, bend

def makeVe_lambda(detVe):
    """- делает Ve в зависимости от detVe
    """

    global btrue
    global leny
    #длина выходного вектора, я не знаю, как её
    x = detVe**(1/leny)
    Ve=np.zeros((leny,leny))
    np.fill_diagonal(Ve,x)
    return Ve

#архитектура этого скрипта такова, что её как будто бы делал наркоман




#INPARAMETERS RANGES
# 'conditions'
# 'plantype'   1 Aprior 0 unif [True,False]
# 'n'  (5,50,5)  5 10 15 20 25 30 35 40 значений
# 'lbinitbtrue'
# 'diapwidth'
# 'assym'
# 'nsiggen'
# 'nsig'
# 'detve'
# 'isbinitgood'




#INPARAMETERS
'conditions'
'plantype'
'n'
'lbinitbtrue'
'diapwidth'
'assym'
'nsiggen'
'nsig'
'detve'
'isbinitgood'



def makeavlst (inl, word):
    return np.average([i[word] for i in inl])


def makeav(minilist):
    keys = list(minilist[0].keys())
    res=dict()
    for key in keys:
        res[key]=makeavlst(minilist, key)
    return res

def fffunctionformatting(res, file):
    printf = lambda x: print (x, file)




#тестить сначала на интервалах, которые единичны
def mainfunc (i_plantype, i_n, i_lbinitbtrue, i_diapwidth, i_assym, i_nsiggen, i_nsig, i_detve,  i_isbinitgood,
              gknuxfunc_lambda,
              makeplan_lambda, makebinit_lambda, makediap_lambda, makeVe_lambda
              ):
    """
        plantype
        n
        lbinitbtrue
        diapwidth
        assym
        nsiggen
        nsig
        detve

        isbinitgood

        Лямбды:
        gknuxfunc_lambda - лямбда оценочной функции, в которой определено всё, кроме того, что задаётся входными параметрами,
        в изменении которых состоит эксперимент

        makeplan_lambda - делает нужный план в зависимости от параметров plantype и n
        makebinit_lambda - делает начальное значение в зависимости от параметров lbinitbtrue
        makediap_lambda - делает диапазон в зависимости от параметров diapwidth, assym
        makeVe_lambda - делает Ve в зависимости от detVe

        btrue фиксировано. Всё остальное определяется, отталкиваясь от него
    """
    res  =  list()
    resfolder= 'Results'
    #проверка доступности папки записи результатов
    with open (resfolder+'/'+'rescsv.csv', 'w')  as file:
        file.write("Commit: ")
        file.write (input("Please a comment, will be written first in results file. Preferably current commit hash: "))
        file.write('\n\n')

    print (locals()) #чтоб не расслабляться при чтении логов

    firstiteration=True

    for plantype in i_plantype:
        for n in i_n:
            for lbinitbtrue in i_lbinitbtrue: #этот цикл идёт по вариантам инит-вектора, который берётся по равномерному распределению из диапазона
                for diapwidth in i_diapwidth: #задаётся в долях от btrue
                    for assym in i_assym: #задаётся в долях от btrue
                        for nsiggen in i_nsiggen: #кандидат на удаление
                            for nsig in i_nsig: #кандидат на удаление
                                for detVe in i_detve: #в случае однооткликовой модели это - единственный элемент оной матрицы
                                    for isbinitgood in i_isbinitgood: #только 1 и 0

                                        bstart, bend = makediap_lambda(diapwidth, assym)
                                        plan = makeplan_lambda(plantype, n, bstart, bend)
                                        Ve = makeVe_lambda (detVe)
                                        bb=makebinit_lambda(bstart, bend)
                                        binit=bb[0]

                                        condition = {'plantype': plantype,
                                                     'n':n,
                                                     'lbinitbtrue':bb[1],
                                                     'diapwidth':diapwidth,
                                                     'assym':assym,
                                                     'nsiggen':nsiggen,
                                                     'nsig':nsig,
                                                     'detve':detVe,
                                                     'isbinitgood':isbinitgood}

                                        dellist = ['log', 'Sklist', 'DispLT', 'DispDif', 'Diflist', 'name', 'Vb', 'VbSigmas']

                                        if detVe<10e-15:  #считаем Ve практически нулевой, разброс - отсутствующим
                                            result = gknuxfunc_lambda (plan, binit, bstart, bend, Ve)
                                            for i in dellist:
                                                del (result[i])

                                        else:
                                            minilist=list()
                                            for i in range (0, 20): #выборка 20
                                                result = gknuxfunc_lambda (plan, binit, bstart, bend, Ve)
                                                for j in dellist:
                                                    del (result[j])
                                                minilist.append(result)
                                            result=makeav(minilist)  #усредняет все показатели в списке
                                        conditionkeys = sorted(condition.keys())
                                        resultkeys = sorted(result.keys())
                                        with open (resfolder+'/'+'rescsv.csv', 'a')  as file:
                                            if firstiteration:
                                                print (','.join(str(x) for x in conditionkeys+resultkeys), file) #вывели заголовок CSV столбцов
                                                firstiteration=False
                                            file.write(','.join(str(condition[x]) for x in conditionkeys))
                                            file.write(','.join(str(result[x]) for x in resultkeys))
                                            file.write('\n')

                                            #данные пишутся постепенно и файл тотчас закрывается, как только они записаны


                                        res.append ({'condition': condition, 'result':result})




    #что мы делаем с res
    #сперва вкатаем в pickle
    with (resfolder+'/'+'respickled0.pkl', 'wb') as f1:
        pickle.dump(res, f1)

    with (resfolder+'/'+'respickledh.pkl', 'wb') as f2:
        pickle.dump(res, f2, -1)

    with (resfolder+'/'+'restext.text', 'wt') as f3:
        fffunctionformatting(res, f3)

    return 0






