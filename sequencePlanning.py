__author__ = 'vasilev_is'
import copy

import numpy as np

import Ofiura.Ofiura_planning as o_p
import Ofiura.Ofiura_Estimation as o_e
import Ofiura.Ofiura_general as o_g


#http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html
def getbSeqPlanUltra (xstart:list, xend:list, N:int, btrue:list, binit:list, c, Ve, jacf, funcf, initplan=None, NSIG=3, smallestdetVb=1e-6, implicit=False, lognorm=False):
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
    :param implicit - неявность функции. Если функция неявна, то в gknu sign=0
    :param lognorm - требуется ли использование логнорм для получения измеренных данных
    :return: k, число итераций, лог
    Для переделывания на реальные измерения btrue=None и функции моделирования переводятся на измерительные
    """
    startplan = initplan if initplan else o_p.makeUniformExpPlan(xstart, xend, N)
    origplan = copy.copy(startplan)
    measdata = o_p.makeMeasAccToPlan_lognorm(funcf, startplan, btrue, c, Ve) if lognorm else o_p.makeMeasAccToPlan(funcf, startplan, btrue, c, Ve)
    log=""
    b=binit
    for numiter in range(100): #ограничитель цикла - если выход произошёл по ограничению, значит, возможна ошибка

        estim=o_e.grandCountGN_UltraX1(funcf, jacf, measdata, b, c, NSIG, implicit=implicit) #получили оценку b binit=b
        b=estim[0]
        Sk=estim[1]
        #Vb=o_ap.countVbForPlan(startplan, b,  c, Ve, jacf, funcf)
        Vb=o_p.countVbForMeasdata(b,  c, Ve, jacf, measdata)
        #посчитали определитель
        detVb=np.linalg.det(Vb)

        print ("iteration: {0}\nb={1}\ndetVb={2}\nSk={3}".format(numiter, b, detVb, Sk))


        #FIXME очень стрёмное условие выхода
        if detVb>0 and detVb<smallestdetVb: #если определитель меньше 10^-6
            return b, numiter, Sk, startplan, origplan, log #то всё вернуть

        #иначе поиск следующей точки плана

        xdot=copy.copy(xstart) #получили начальную точку начальное значение - ровно пополам диапазона
        for i in range(len(xstart)): #присвоили ей значение
            xdot[i]=xstart[i]+(xend[i]-xstart[i])/2

        #объектная функция
        function = lambda x: np.linalg.det(o_p.countVbForPlan(o_g.appendToList(startplan, x),b,c,Ve,jacf, funcf))
        #FIXME здесь важный момент - под func обычно понимается бездисперсионная функция, что в данном случае неверно.
        #function и measure есть разные функции - первая даёт идеальный результат, вторая - с дисперсией
        #создать функцию, которая будет возвращать полученные от func данные, налагая дисперсию.

        #каждый раз будет пытаться добавить в план точку и вернуть определитель с добавленной точкой

        xdot=o_g.doublesearch (xstart, xend, xdot, function) #оптимизировали значение точки

        startplan.append(xdot)
        measdata.append({'x':xdot, 'y':funcf(xdot,b,c)})

    #окончание этого цикла "естественным путём" говорит о том, что превышено максимальное число итераций
    return b, 100, Sk, startplan, origplan, log+"ERROR: maximum number of iterations archieved" #то всё вернуть


def test():
    """
    Тестирует последовательное планирование
    :return:
    """
    xstart=[1, 100]
    xend=[20,200]

    N=10
    c={"a":1000}
    funcf=lambda x,b,c: np.array ( [ b[0]+b[1]*x[0]+b[2]*x[1]+b[3]*x[0]*x[1]+b[4]*x[0]*x[0]+b[5]*x[1]*x[1],   b[6]+b[7]*x[0]+b[8]*x[1]+b[9]*x[0]*x[1]+b[10]*x[0]*x[0]+b[11]*x[1]*x[1] ] )
    jacf = lambda x,b,c,y: np.matrix([ [1, x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1], 0, 0, 0, 0, 0, 0],
                                      [0,0,0,0,0,0,1,x[0], x[1], x[0]*x[1], x[0]*x[0], x[1]*x[1]] ]).T

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

    print ("b, numiter, Sk, plan, origplan, log")
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

