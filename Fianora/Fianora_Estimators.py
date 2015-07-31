__author__ = 'vasilev_is'

import  Fianora.Fianora_Models
import numpy as np
import platform
import copy
import math
from prettytable import PrettyTable



class Options():
    """
    Estimation process options
    """

    def __getattr__(self, item):
        return None

class AbstractEstimator():

    def init_parameters(self, ec, model):
        """
        Эта функция должна вызываться у объекта всегда перед оценкой
        :param bstart:
        :param bend:
        :param binit:
        :param xstart:
        :param xend:
        :param model:
        :param measdata:
        :param Ve:
        :return:
        """

        self.bstart = ec.bstart
        self.bend = ec.bend
        self.binit = ec.binit
        self.xstart = ec.xstart
        self.xend = ec.xend
        self.model = model
        self.Ve = ec.Ve

    def estimate (self, measdata, options = None):
        """
        Шаблонный метод, идея в том, что сначала выполняется оценка, потом определяются показатели адекватности
        в итоге получается некоторый стандартный словарь
        :return:
        """
        z = self.estimate_method(measdata, options).copy()
        z.update(self.adequacy(z, measdata))
        z['name'] = self.model.name

        return z

    def estimate_method(self, measdata, options):
        """
        Estimation method, like gknux
        :return:
        """
        raise NotImplementedError ("Pleas implement me")

    def adequacy(self, est, measdata):
        """
        Adequacy determiner
        Спецификация:
        Тип данных - словарь
        Поля и описание -
         'b'            оцененный вектор коэффициентов
         'numiter'      число итераций
         'log'          лог
         'Sklist'       список значений объектной функции
         'Sk'           последнее значение объектной функции

         'Vb'           ковариационная матрица оценённого вектора b
         'VbSigmas'     сигмы из ковариационной матрицы оценённого вектора b (корни квадратные из диагонали)

         'AvLogTruth'   среднее по точкам эксперимента значение логарифма правдоподобия
         'DispLT'       дисперсия логарифма правдоподобия
         'SigmaLT'      среднеквадратическое отклонение логарифма правдоподобия
         'AvDif'        средний остатков
         'DispDif'      дисперсия остатков
         'SigmaDif'     СКВ остатков
         'Diflist'      список остатков
         'name'         название метода
        """

        rs={}

        b = est['b']
        Vb = rs['Vb'] = self.countVbForMeasdata(b, measdata)   # np.linalg.inv(G)

        sigmas=list()
        for i in range (Vb.shape[0]):
            sigmas.append(math.sqrt(Vb[i][i]))
        rs['VbSigmas'] = sigmas

        logTruthnessStuff = self.logTruthness(b, measdata) #  Average, Disp, math.sqrt(Disp)
        diflistStuff =  self.averageDif(b, measdata) # np.average(diflistn), np.var(diflistn), math.sqrt(np.var(diflistn))
        names=['AvLogTruth','DispLT', 'SigmaLT', 'AvDif', 'DispDif', 'SigmaDif', 'Diflist']
        values =  list(logTruthnessStuff) + list(diflistStuff)
        rs.update(dict(zip (names, values)))


        return rs

    def countVbForMeasdata(self, b:list, measdata):
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

        bstart = self.bstart
        bend = self.bend
        binit = self.binit
        xstart = self.xstart
        xend = self.xend
        model = self.model

        Ve = self.Ve
        func = self.model.funcf
        jac = self.model.jacf

        G=np.zeros((len(b),len(b))) #матрица G


        for point in measdata:
            jj=jac(point['x'], b)
            #G+=jj*np.linalg.inv(Ve)*jj.T

            #print  (jj.T, np.linalg.inv(Ve), jj)

            G+=np.dot ( np.dot(jj.T, np.linalg.inv(Ve)), jj)

            #G+=np.dot(jj.T, jj)
        try:
            return np.linalg.inv(G)
        except BaseException as e:
            print('Fatal error in countVbForMeasdata: ',e)
            print('b vector=',b)
            print('current point=',point)
            print('G=',G)
            exit(0)

    def logTruthness (self, b, measdata):
        """
        Считает значение объектной функции (!) это не логарифм функции правдоподобия!
        :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
        :param b: оценка вектора коэффициентов
        :param Ve: ковариационная матрица ошибок экспериментальных данных
        :param func: callable функция x,b,c, возвращает значение y
        :param c: словарь дополнительных переменных
        :return: среднее значение логарифма функции правдоподобия, дисперсию по выборке экспериментальных данных, стандартное отклонение
        """

        bstart = self.bstart
        bend = self.bend
        binit = self.binit
        xstart = self.xstart
        xend = self.xend
        model = self.model

        Ve = self.Ve
        func = self.model.funcf

        S=list()

        for i in range(len(measdata)):
            measpoint = measdata[i]
            dif=np.array(measpoint['y'])-np.array(func(measpoint['x'],b))

            S.append(np.dot(np.dot(dif.T, np.linalg.inv(Ve)), dif))

        K=Ve.shape[0] #число откликов
        M=len(b) #число коэффициентов
        N=len(measdata)
        shift=K*N/(K*N-M)

        Snp=list(map(np.float64, S))

        Average=np.average(Snp)*shift
        Disp = np.var(S)*shift*shift

        return Average, Disp, math.sqrt(Disp)

    def averageDif(self, b, measdata):
        """
        :param measdata: список словарей экспериментальных данных [{'x': [] 'y':[])},{'x': [] 'y':[])}]
        :param b: вектор коэфф
        :param Ve:  Ve ковар. матрица измеренных данных
        :param funcf callable функция, параметры по формату x,b,c
        :param c словарь дополнительных постоянных
        :return: среднее, дисперсия, стандартное отклонение
        """


        bstart = self.bstart
        bend = self.bend
        binit = self.binit
        xstart = self.xstart
        xend = self.xend
        model = self.model

        Ve = self.Ve
        func = self.model.funcf


        diflist=list()
        diflistNoAbs=list()


        for measpoint in measdata:
            dif=np.array(measpoint['y'])-func(measpoint['x'],b)

            diflist.append(np.abs(dif ))
            diflistNoAbs.append(dif)

            diflistn=list(map(np.float64, diflist))
            diflistna=list(map(np.float64, diflistNoAbs))


        return np.average(diflistn), np.var(diflistn), math.sqrt(np.var(diflistn)), diflistna


class NGEstimator(AbstractEstimator):
    """
    Оценщик, использующий метод Гаусса-Ньютона с ограничениями
    """

    def makeSkInit (self, measdata):

        funcf = self.model.funcf
        binit = self.binit
        Sk=None
        for point in measdata:
            fxbc=funcf(point['x'],binit)
            if fxbc is None:
                print ("ERROR: makeSkInit: funcf returned None")
                return None
            dif=point['y']-fxbc
            if Sk is None:
                Sk=np.dot(dif.T,dif)
            else:
                Sk+=np.dot(dif.T,dif)
            return Sk

    def islinux (self):
        if 'Windows' in platform.system():
            return False
        return True

    def makeAinit (self,  measdata, isBinitGood=True):
        """
        расчёт ведётся по правилу "binit есть хорошее предположение в рамках области"
        :param bstart:
        :param bend:
        :param Skbinit:
        :param binit:
        :return:
        """
        bstart = self.bstart
        binit = self.binit
        bend = self.bend
        Skbinit = self.makeSkInit(measdata)

        A=np.zeros ( (len(binit), 2 ))  #,  dtype=np.float96 if islinux() else np.float64

        for i in range (len(binit)):
            #A[i][0]=A[i][1]=.001*(bend[i]-bstart[i])*Skbinit  #универсально

            if isBinitGood:
                A[i][0] = 0.001*(binit[i]-bstart[i])*Skbinit
                A[i][1] = 0.001*(bend[i]-binit[i])*Skbinit
            else:
                A[i][0] = A[i][1] = 0.001*(bend[i]-bstart[i])*Skbinit  #универсально
        return A

    def countSklims(self, A, b):
        bstart = self.bstart
        bend = self.bend
        partone=parttwo=.0
        for i in range (len(b)):
            partone+=A[i][0]/(b[i]-bstart[i]) #UNSAFE: если вдруг, ну вдруг b=bstart или bend, то будет <strike> треш </strike> креш с делением на ноль.
            parttwo+=A[i][1]/(bend[i]-b[i])
        return partone+parttwo

    def countN (self, A, b):
        bstart = self.bstart
        bend = self.bend

        N=np.zeros ((len(b),len(b)))
        for j in range (len(b)):
            partone=2*A[j][0]/(b[j]-bstart[j])**3
            parttwo=2*A[j][1]/(bend[j]-b[j])**3
            N[j][j]+=parttwo+partone #так как матрица нулевая
        return N

    def estimate_method(self, measdata, options):
        binit = self.binit

        NSIG = options.NSIG if options.NSIG else 50
        NSIGGENERAL = options.NSIGGENERAL if options.NSIGGENERAL else 50
        verbose = options.verbose if options.verbose else False
        verbose_wrapper = options.verbose_wrapper if options.verbose_wrapper else False
        isBinitGood = options.isBinitGood if options.isBinitGood else False



        if isinstance(self.model, Fianora.Fianora_Models.ImplicitModel): #если модель - это подкласс неявной модели, то...
            implicit = True
        else:
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

        A=self.makeAinit(measdata, isBinitGood)

        log=''

        if verbose_wrapper:
            print ('==grandCountGN_UltraX1_Limited_wrapper is launched==\n\n')

        for numiter in range (maxiter):
            bpriv=copy.copy(b)

            gknux=self.grandCountGN_UltraX1_Limited (A, measdata, NSIG, implicit, verbose) #посчитали b


            if gknux is None:
                print ("grandCountGN_UltraX1_Limited_wrapper crashed on some iteration")
                return None
            gknuxlist.append(gknux)
            if verbose_wrapper:
                print ('Iteration \n',numiter,'\n' ,gknux)

            b=gknux['b']


            if not gknux['log']=='':
                #log+="On gknux iteration "+numiter+": "+ gknux[2]+"\n"
                log+="On gknux iteration {0}: {1}\n".format (numiter, gknux[2])


            for j in range (len(binit)): #уменьшили в два раза
                A[j][0]*=0.5
                A[j][1]*=0.5

            condition=False
            for i in range (len(b)):
                if math.fabs ((b[i  ]-bpriv[i])/bpriv[i]) > math.pow(10,-1*NSIGGENERAL):
                    condition=True
            if not condition:
                break

            #мол если хоть один компонент вектора b значимо изменился, тогда продолжать. Иначе программа дойдёт до break и цикл прекратится

        if verbose_wrapper:
            print ('grandCountGN_UltraX1_Limited_wrapper iterations number:', numiter)
            print (gknux, log)

        gknux['log'] = log
        return gknux

    def grandCountGN_UltraX1_Limited (self, A, measdata,  _NSIG=3, implicit=False, verbose=False):
        """
        Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом с ограничениями, заданными диапазоном
        В стандартный поток вывода выводит отладочную информацию по каждой итерации


        :param NSIG=3 точность (кол-во знаков после запятой)
        :param implicit True если функция - неявная, иначе false
        :param verbose Если True, то подробно принтить результаты итераций
        :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
        """
        funcf = self.model.funcf
        jacf = self.model.jacf

        binit = self.binit

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
                jac=jacf(point['x'],b)

                if jac is None:
                    return None
                G+=np.dot(jac.T,jac)

                try:
                    dif=np.array(point['y'])-np.array(funcf(point['x'],b))
                except BaseException as e:
                    if verbose:
                        print('grandCountGN_UltraX1_limited: As funcf returned None, method  stops:', e)
                        print (point, b, sep='\n')
                    return None


                if B5 is None:
                    B5=np.dot(dif, jac)
                else:
                    B5+=np.dot(dif,jac)
                Sk+=np.dot(dif.T,dif)


            N=self.countN(A,b)
            Sklims=self.countSklims(A,b)

            #G=G-N if implicit else G+N #добавляем градиент от штрафных функций
            G=G+N
            Sk+=Sklims #добавиляем объектную функцию от штрафных функций

            #print(np.linalg.inv(G), B5[:,0])
            #костыль для диодной задачи
            if hasattr(B5, 'A1'):
                B5=B5.A1
            try:
                deltab=np.dot(np.linalg.inv(G), B5)
            except BaseException as e:
                if verbose:
                    print('Error in G:', e)
                    print('G=',G)
                return None


            #mu counting
            mu=4
            cond2=True
            it=0
            Skmu=0
            while (cond2):
                Skmu=self.countSklims(A,b) #добавляем объектную функцию от штрафных функций
                mu/=2
                for point in measdata:
                    try:
                        dif=np.array(point['y'])-np.array(funcf(point['x'],b+deltab*mu)) if sign else np.array(point['y'])-np.array(funcf(point['x'],b-deltab*mu))
                    except:
                        continue

                    Skmu+=np.dot(dif.T, dif)

                it+=1


                if (it>200):
                    log+="Mu counting: break due to max number of iteration exceed"
                    break
                cond2=Skmu>Sk

            b=b+deltab*mu if sign else b-deltab*mu

            Sk=Skmu

            Sklist.append(Sk)
            if verbose:
                print ("Sk:",Sk)
                print ("Iteration {0} mu={1} delta={2} deltamu={3} resb={4}".format(numiter, mu, deltab, deltab*mu, b))

            numiter+=1

            condition=False
            for i in range (len(b)):
                if math.fabs ((b[i]-bpriv[i])/bpriv[i])>math.pow(10,-1*_NSIG):
                    condition=True


            if numiter>3000: #max number of iterations
                log+="GKNUX1: Break due to max number of iteration exceed"
                break

        if len(Sklist)>100:
            ll=[Sklist[i] for i in range(0,len(Sklist),10)]
            Sklist = ll

        return {'b':b, 'numiter':numiter, 'log':log, 'Sklist':Sklist, 'Sk':Sk}


class AbstractEstimatorDecorator(AbstractEstimator):
    def __init__(self, component:AbstractEstimator):
        self.component = component

    def estimate(self, measdata, options):
        #поведение по умолчанию
        return self.component.estimate(options)


class ConsoleEstimatorDecorator(AbstractEstimatorDecorator):
    def estimate(self, measdata, options):
        """
        Усложнённый декоратор: во-первых, параметры пробрасываются от конца к началу цепочки,
        во-вторых, результат пробрасываться от начала к концу цепочки
        :param options:
        :return:
        """
        est = self.component.estimate(measdata, options)

        tg=copy.copy(est) #для вывода в таблицу

        del tg['log']
        del tg['Sklist']
        del tg['Diflist']
        del tg['name']
        del tg['Vb']

        klist = sorted(list(tg.keys()))
        vallist = [tg[i] for i in klist]



        num=2 #здесь некоторая магия. Дефолтно - 4, текущие -2, и, вроде, показывает все значения
        klists=[klist[i::num] for i in range(num)]
        vallists=[vallist[i::num] for i in range(num)]


        t=PrettyTable (klists[0])
        t.add_row(vallists[0])

        t2=PrettyTable (klists[1])
        t2.add_row(vallists[1])


        print ("\n Results of method {0} \n".format(est['name']))
        print (t)
        print (t2)
        print ("Vb:\n {0}  ".format (est['Vb']))
        print ("\nLog messages: {0} \n ".format (est['log']))



        return est




