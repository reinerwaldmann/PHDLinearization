__author__ = 'vasilev_is'

import pickle
import hashlib
import time
import datetime

import numpy as np

import Ofiura.Ofiura_general as o_g
import Ofiura.Ofiura_ApriorPlanning  as o_ap
import Ofiura.Ofiura_EstimationLimited  as o_el
import Cases.Dependency.Diode_In_Limited_Dependency as cdld
import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Qualitat as o_q




#конкретно к модели диода относящиеся вещи

btrue = [1.238e-14, 1.8, 100]
c={}
leny=1 #длина выходного вектора
xstart, xend =   [0.001], [1.4]

funcf = cdld.solver_Diode_In_Limited
jacf = cdld.jac_Diode_In_Limited


usercomment=''

#print (np.linalg.norm(btrue))



def  gknuxfunc_lambda (plan, binit, bstart, bend, Ve):
    """

    :param plan:
    :param binit:
    :param bstart:
    :param bend:
    :param Ve:
    :return:
    """
    global btrue, c, xstart, xend, funcf, jacf
    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, plan, btrue, c,Ve )
    gknuxlim = o_el.grandCountGN_UltraX1_Limited_wrapper(funcf,jacf,measdata,binit,bstart,bend, c, implicit=True, verbose=False, verbose_wrapper=False )
    try:
        gknuxlim2=o_q.convertToQualitatStandart (gknuxlim, funcf, jacf,  measdata, c, Ve, name='Limited Count Aprior Dependency')
    except:
        return None
    return gknuxlim2


def makeplan_lambda (plantype, n, bstart, bend, Ve):
    global xstart, xend, btrue,c
    #- делает нужный план в зависимости от параметров plantype и n
    #Надо создать пачку априорных планов - 5 10 15 20 25 30 35 40 значений
    hash= hashlib.md5((n,bstart,bend,Ve).__str__().encode('utf-8')).hexdigest()

    #hashlib.md5("whatever your string is".encode('utf-8')).

    filename = usercomment+'_'+str(hash)
    if plantype: #если априорный
        return o_ap.makePlanCached (xstart, xend, n, bstart, bend, c, Ve, jacf, funcf, Ntries=4, verbose=False, foldername='DependencyPlans', cachname=filename)
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

    global usercomment
    res  =  list()
    resfolder= 'Results'
    #проверка доступности папки записи результатов
    with open (resfolder+'/'+'rescsv.csv', 'wt')  as file:
        file.write("Commit: ")
        usercomment=input("Please a comment, will be written first in results file. Preferably current commit hash. Mind that cached plans will be prefixed with this line: ")
        file.write (usercomment)
        file.write('\n\n')

#    print (locals()) #чтоб не расслабляться при чтении логов

    firstiteration=True
    iternum=0

    #for plantype in i_plantype:
    for plantype in [1]:
        for n in i_n:
            for lbinitbtrue in i_lbinitbtrue: #этот цикл идёт по вариантам инит-вектора, который берётся по равномерному распределению из диапазона
                for diapwidth in i_diapwidth: #задаётся в долях от btrue
                    for assym in i_assym: #задаётся в долях от btrue
                        for nsiggen in i_nsiggen: #кандидат на удаление
                            for nsig in i_nsig: #кандидат на удаление
                                for detVe in i_detve: #в случае однооткликовой модели это - единственный элемент оной матрицы
                                    for isbinitgood in i_isbinitgood: #только 1 и 0

                                        prevtime = time.time()
                                        print ('iteration ', iternum)



                                        bstart, bend = makediap_lambda(diapwidth, assym)
                                        Ve = makeVe_lambda (detVe)
                                        plan = makeplan_lambda(plantype, n, bstart, bend,Ve)
                                        bb=makebinit_lambda(bstart, bend)
                                        binit=bb[0]

                                        print ('mainfunc: iteration parameters:')
                                        print (iternum, binit, bstart, bend)

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
                                            if result is not None:
                                                for i in dellist:
                                                    del (result[i])
                                        else:
                                            minilist=list()
                                            for i in range (0, 20): #выборка 20
                                                result = gknuxfunc_lambda (plan, binit, bstart, bend, Ve)

                                                if result is not None:
                                                    for j in dellist:
                                                        del (result[j])
                                                    minilist.append(result)

                                            if (len(minilist)>0):
                                                result=makeav(minilist)  #усредняет все показатели в списке
                                            else:
                                                result=None

                                        conditionkeys = sorted(condition.keys())
                                        resultkeys = sorted(result.keys())
                                        with open (resfolder+'/'+'rescsv.csv', 'a')  as file:
                                            if iternum==0:
                                                file.write(','.join(str(x) for x in ['iteration num',]+conditionkeys+resultkeys)) #вывели заголовок CSV столбцов
                                                file.write('\n')

                                            file.write(str(iternum).join(','))
                                            file.write(','.join(str(condition[x]) for x in conditionkeys))


                                            if  result is not None:
                                                file.write(','.join(str(result[x]) for x in resultkeys))
                                            else:
                                                file.write(','.join('None' for x in resultkeys))
                                            file.write('\n')

                                            #данные пишутся постепенно и файл тотчас закрывается, как только они записаны
                                        res.append ({'condition': condition, 'result':result})

                                        st = datetime.datetime.fromtimestamp(time.time()-prevtime).strftime('%H:%M:%S')
                                        print ('Iteration number {0} finished in time {1}'.format(iternum, st  ))
                                        iternum+=1


    #что мы делаем с res
    #сперва вкатаем в pickle
    with (resfolder+'/'+'respickled0.pkl', 'wb') as f1:
        pickle.dump(res, f1)

    with (resfolder+'/'+'respickledh.pkl', 'wb') as f2:
        pickle.dump(res, f2, -1)

    return 0



def test ():
    """
    Вся необходимая информация для запуска размещается либо в глобальных переменных, на которые ссылаются функции напрямую, либо в Diode_In_LimitedDependency,
    на которые идут ссылки из тех же глобальных переменных
    Mainfunc launcher
    :return: nothing
    """
    #RANGES:
    #Здесь задаются диапазоны изменения входных параметров эксперимента


    i_plantype=[1] #мы будем  работать только с априорными планами, их эффективность можно доказать и с помощью более узких экспериментов
    #но эксперименты по сравнению нужны!
    #априорный план может сработать плохо при очень широких границах, возможно, нужен-таки равномерный план
    i_n=[10,20,30,40] #число точек в плане (априорном, внимание!)
    i_lbinitbtrue = range(10) #десять попыток uniform-выбора
    i_diapwidth = np.arange(0.10, 0.4, 0.05) #ширина диапазона от десяти процентов до 40 процентов с шагом в 5 процентов //6
    i_assym = np.arange(0.001, 0.04, 0.01) #ассиметрия //
    i_nsiggen = (20,)
    i_nsig = (20,)
    i_detve = (10e-2, 10e-3, 10e-4, 10e-7, 10e-10, 10e-20) #последнее значение отключает внедрение дисперсии, то есть y становится неслучайным
    i_isbinitgood =[0,1] #включать или нет "хорошесть" начального приближения A

    #вывод параметров поставленной задачи - длина всех последовательностей
    tl=1
    for i, val in locals().items():
        if i.startswith('i_'):
            print (i, len(val))
            tl*=len(val)
    print ('Total task length is {0} iterations'.format(tl))


    #LAMBDAS
    mainfunc (i_plantype, i_n, i_lbinitbtrue, i_diapwidth, i_assym, i_nsiggen, i_nsig, i_detve,  i_isbinitgood,
              gknuxfunc_lambda,
              makeplan_lambda, makebinit_lambda, makediap_lambda, makeVe_lambda
              )


test ()




#список более узких экспериментов:
#2 кривые: зависимость качества от количества точек в плане, одна кривая для равномерного плана, другая - для априорного