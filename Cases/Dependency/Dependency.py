__author__ = 'vasilev_is'



import pickle

import numpy as np
import Ofiura.Ofiura_general as o_g


#конкретно к модели диода относящиеся вещи

btrue = [1.238e-14, 1.8, 100]
#print (np.linalg.norm(btrue))


def makeplan_l (plantype, n ):
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



def makediap_lambda (diapwidth, assym):
    #- делает диапазон в зависимости от параметров diapwidth, assym
    global btrue

    btrue1=np.array(btrue)

    bcenter = np.array(btrue) * (1+np.linalg.norm(btrue)/assym)

    print (np.linalg.norm(btrue1)/assym)

    bstart, bend = bcenter* (1-np.linalg.norm(bcenter)/(2*diapwidth)), bcenter* (1+np.linalg.norm(bcenter)/(2*diapwidth))
    return bstart, bend

print (makediap_lambda (0.1, 0.2))



def makeVe_lambda(detVe):
    #- делает Ve в зависимости от detVe
    global btrue



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
    with open (resfolder+'/'+'test.txt', 'wt')  as file:
        file.write('test')

    print (locals()) #чтоб не расслабляться при чтении логов

    for plantype in i_plantype:
        for n in i_n:
            for lbinitbtrue in i_lbinitbtrue:
                for diapwidth in i_diapwidth:
                    for assym in i_assym:
                        for nsiggen in i_nsiggen:
                            for nsig in i_nsig:
                                for detVe in i_detve:
                                    for isbinitgood in i_isbinitgood:

                                        binit=bb[0]
                                        plan = makeplan_lambda(plantype, n)
                                        bstart, bend = makediap_lambda(diapwidth, assym)
                                        Ve = makeVe_lambda (detVe)
                                        bb=makebinit_lambda(bstart, bend)

                                        condition = {'plantype': plantype,
                                                     'n':n,
                                                     'lbinitbtrue':bb[1],
                                                     'diapwidth':diapwidth,
                                                     'assym':assym,
                                                     'nsiggen':nsiggen,
                                                     'nsig':nsig,
                                                     'detve':detVe,
                                                     'isbinitgood':isbinitgood}


                                        if detVe<10e-15:  #считаем Ve практически нулевой, разброс - отсутствующим
                                            result = gknuxfunc_lambda (plan, binit, bstart, bend, Ve)
                                            del  (result['log'])
                                            del  (result['Sklist'])
                                            del  (result['DispLT'])
                                            del  (result['DispDif'])
                                            del  (result['Diflist'])
                                            del  (result['name'])
                                            del  (result['Vb'])
                                            del  (result['VbSigmas'])
                                        else:
                                            minilist=list()
                                            for i in range (0, 20): #выборка 20
                                                result = gknuxfunc_lambda (plan, binit, bstart, bend, Ve)
                                                del  (result['log'])
                                                del  (result['Sklist'])
                                                del  (result['DispLT'])
                                                del  (result['DispDif'])
                                                del  (result['Diflist'])
                                                del  (result['name'])
                                                del  (result['Vb'])
                                                del  (result['VbSigmas'])
                                                minilist.append(result)
                                            result=makeav(minilist)

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






