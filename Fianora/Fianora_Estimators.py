"""
Описание выходной структуры данных любого оценщика или его декоратора
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
         'GammaDelta'   радиус эллипса в точке оптимума

"""
__author__ = 'vasilev_is'

import  Fianora.Fianora_Models
import numpy as np
import platform
import copy
import math
from prettytable import PrettyTable
import Fianora.Fianora_static_functions as f_sf




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
        #  функция adequacy здесь - какбы функция - декоратор, внутренняя
        #  есть структура данных,est, которая передаётся функции(ям), которые дополняют её на основе своих
        #  сведений и сведений, определённых в ней, некоторыми данными

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
         'GammaDelta'   полширины эллипсоида рассеяния по измерению данного элемента вектора b
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
        rs['GammaDelta'] = self.countGammaDelta(b, measdata, est['Sk'])


        return rs

    def countGammaDelta (self, b:list, measdata, Sk):
        """
        P174 Bard.

        b_i*-GammaDelta<=b_i<=b_i*+GammaDelta

        GammaDelta = (2e/A[i][i])^.5
        where
        A = H*
        H* - Гессиан функции правдоподобия
        Считаем, что H* = sum(jj(measdata(i,b*).T*jj(measdata(i,b*))
        то есть сумма по точкам
        :param b:
        :param measdata:
        :return: вектор гамма-дельта значений,
        равных полширине эллипсоида рассеяния
        """
        jac = self.model.jacf

        H=np.zeros((len(b),len(b))) #матрица G

        for point in measdata:
            jj=jac(point['x'], b, point['y'])
            H+=np.dot(jj.T, jj)



        H*=2
        e=Sk
        return [math.sqrt(2*e/H[i][i]) for i in range (len(b)) ]







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
            jj=jac(point['x'], b, point['y'])
            G+=np.dot ( np.dot(jj.T, np.linalg.inv(Ve)), jj)
            #G+=np.dot ( jj.T, jj)



        try:
            return np.absolute(np.linalg.inv(G))
            #return np.absolute(G)
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


        diflist = [np.array(measpoint['y'])-func(measpoint['x'],b) for measpoint in measdata]
        if len(diflist[0])>0:
            #  если риальни этот список являет собой список векторов, то есть модель многооткликовая
            Tina_LaFee_list = [np.average([x[i] for x in diflist if x[i] < 1e10]) for i in range (len(diflist[0]))] #  average diflist
            varlist = [np.var([x[i] for x in diflist if x[i] < 1e10]) for i in range (len(diflist[0]))] #  var diflist
            sigmaa = list( map (math.sqrt, varlist))
        else:
            Tina_LaFee_list = np.average(diflist)
            varlist = np.var(diflist)
            sigmaa = math.sqrt(varlist)

        return Tina_LaFee_list, varlist, sigmaa, diflist

class NGEEstimatorUnlimited (AbstractEstimator):

    def estimate_method(self, measdata, options):
        """
        Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом
        В стандартный поток вывода выводит отладочную информацию по каждой итерации
        :param NSIG=3 точность (кол-во знаков после запятой)
        :param implicit True если функция - неявная, иначе false

        :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
        """
        print ('Estimated Unlimited')

        binit = self.binit
        NSIG = options.NSIG if options.NSIG else 8
        #NSIG определяется рабочей точностью вычислителя. Матрицы мы считаем с точностью 9 значащих, возьмём один в запас
        verbose = options.verbose if options.verbose else False
        if isinstance(self.model, Fianora.Fianora_Models.ImplicitModel): #если модель - это подкласс неявной модели, то...
            implicit = True
        else:
            implicit = False


        # маразм! но как сделать лучше: объект-контейнер может быть None, и каждый его член тоже может быть None
        lst_data_abt_iteration = None
        if options is not None:
            lst_data_abt_iteration = options.lst_data_abt_iteration if options.lst_data_abt_iteration else None


        funcf = self.model.funcf
        jacf = self.model.jacf
        binit = self.binit

        if lst_data_abt_iteration:
            lst_data_abt_iteration.accept(numiter=0, b=binit)

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
                Skmu=0
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
                if math.fabs ((b[i]-bpriv[i])/bpriv[i]) > math.pow(10,-1*NSIG):
                    condition=True

            if lst_data_abt_iteration:
                lst_data_abt_iteration.accept(numiter=numiter, b=b)

            if numiter>50: #max number of iterations
                log+="GKNUX1: Break due to max number of iteration exceed"
                break

        if len(Sklist)>100:
            ll=[Sklist[i] for i in range(0,len(Sklist),10)]
            Sklist = ll

        return {'b':b, 'numiter':numiter, 'log':log, 'Sklist':Sklist, 'Sk':Sk}





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

        NSIG = options.NSIG if options.NSIG else 8
        NSIGGENERAL = options.NSIGGENERAL if options.NSIGGENERAL else 8

        #NSIG определяется рабочей точностью вычислителя. Матрицы мы считаем с точностью 9 значащих, возьмём один в запас

        verbose = options.verbose if options.verbose else False
        verbose_wrapper = options.verbose_wrapper if options.verbose_wrapper else False
        isBinitGood = options.isBinitGood if options.isBinitGood else False

        lst_data_abt_iteration = options.list_data_abt_iteration if options.list_data_abt_iteration else None
        # объект, имеющий свойство Accept - туда отправляют нечто, он это записывает в список


        if isinstance(self.model, Fianora.Fianora_Models.ImplicitModel): #если модель - это подкласс неявной модели, то...
            implicit = True
        else:
            implicit = False

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

        maxiter=50
        b,bpriv=binit,binit
        gknux=None
        gknuxlist=list()

        A=self.makeAinit(measdata, isBinitGood)

        log=''

        if verbose_wrapper:
            print ('==grandCountGN_UltraX1_Limited_wrapper is launched==\n\n')

        maxiter = 1
        for numiter in range (maxiter):
            bpriv=copy.copy(b)

            if verbose_wrapper:
                with f_sf.Profiler() as p:
                   gknux=self.grandCountGN_UltraX1_Limited (A, measdata, NSIG, implicit, verbose, options) #посчитали b
            else:
                gknux=self.grandCountGN_UltraX1_Limited (A, measdata, NSIG, implicit, verbose, options) #посчитали b


            if gknux is None:
                print ("grandCountGN_UltraX1_Limited_wrapper crashed on some iteration")
                return None
            gknuxlist.append(gknux)
            if verbose_wrapper:
                print ('Iteration \n',numiter,'\n' ,gknux)

            b=gknux['b']


            if not gknux['log']=='':
                #log+="On gknux iteration "+numiter+": "+ gknux[2]+"\n"
                log+="On gknux iteration {0}: {1}\n".format (numiter, gknux['log'])


            for j in range (len(binit)): #уменьшили в два раза
                A[j][0]*=0.5
                A[j][1]*=0.5

            condition=False
            for i in range (len(b)):
                if math.fabs ((b[i]-bpriv[i])/bpriv[i]) > math.pow(10,-1*NSIGGENERAL):
                    condition=True
            if not condition:
                break

            #мол если хоть один компонент вектора b значимо изменился, тогда продолжать. Иначе программа дойдёт до break и цикл прекратится

        if verbose_wrapper:
            print ('grandCountGN_UltraX1_Limited_wrapper iterations number:', numiter)
            print (gknux, log)

        gknux['log'] = log
        return gknux

    def grandCountGN_UltraX1_Limited (self, A, measdata,  _NSIG, implicit=False, verbose=False, options=None):
        """
        Производит оценку коэффициентов по методу Гаусса-Ньютона с переменным шагом с ограничениями, заданными диапазоном
        В стандартный поток вывода выводит отладочную информацию по каждой итерации


        :param NSIG=3 точность (кол-во знаков после запятой)
        :param implicit True если функция - неявная, иначе false
        :param verbose Если True, то подробно принтить результаты итераций
        :returns b, numiter, log - вектор оценки коэффициентов, число итераций, сообщения
        """



        # маразм! но как сделать лучше: объект-контейнер может быть None, и каждый его член тоже может быть None
        lst_data_abt_iteration = None
        if options is not None:
            lst_data_abt_iteration = options.lst_data_abt_iteration if options.lst_data_abt_iteration else None




        funcf = self.model.funcf
        jacf = self.model.jacf

        binit = self.binit

        if lst_data_abt_iteration:
            lst_data_abt_iteration.accept(numiter=0, b=binit)



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
                if math.fabs ((b[i]-bpriv[i])/bpriv[i]) > math.pow(10,-1*_NSIG):
                    condition=True

            if lst_data_abt_iteration:
                lst_data_abt_iteration.accept(numiter=numiter, b=b)

            if numiter>50: #max number of iterations
                log+="GKNUX1: Break due to max number of iteration exceed"
                break

        if len(Sklist)>100:
            ll=[Sklist[i] for i in range(0,len(Sklist),10)]
            Sklist = ll

        return {'b':b, 'numiter':numiter, 'log':log, 'Sklist':Sklist, 'Sk':Sk}


class AbstractEstimatorDecorator(AbstractEstimator):
    def __init__(self, component:AbstractEstimator):
        self.component = component

    def init_parameters(self, ec, model):
        self.component.init_parameters(ec, model)

    def estimate(self, measdata, options):
        #поведение по умолчанию
        return self.component.estimate(options)


class ConsoleEstimatorDecorator(AbstractEstimatorDecorator):
    def estimate(self, measdata, options):
        """
        Усложнённый декоратор: во-первых, параметры пробрасываются от конца к началу цепочки,
        во-вторых, результат пробрасываться от начала к концу цепочки

        Конкретно в этом, консольном декораторе делается так - выкидываются лишние поля и остальное выводится в таблице
        лишние поля (см ниже, где del) выводятся потом особо

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

        #[print (f) for f in est['Diflist']]
        #[print (f) for f in measdata]




        return est


class GraphPackEstimatorDecorator (AbstractEstimatorDecorator):
    """
    Декоратор, который либо сохраняет рисунки, либо выдаёт последовательно на экран
    рисунки такие:
        два графика (риальни, моделированни)
		диаграмма остатков
		график падения объектной функции
    """

    def __init__(self, component, foldername, appendix):
        """

        :param component:
        :param foldername: папка, куда сохранять, если None, то показывать
        :param appendix: добавлять к именам файлов
        :return:
        """
        AbstractEstimatorDecorator.__init__(self, component)
        self.foldername = foldername
        self.appendix = appendix

    def estimate(self, measdata, options):
        """
        Описание выходной структуры данных любого оценщика или его декоратора
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
        import matplotlib.pyplot as plt

        est = self.component.estimate(measdata, options)

        try:
            for i in range (len( est['Diflist'][0])):

                part_diflist = [f[i] for f in est['Diflist']]

                part_diflist = [x for x in part_diflist if x<1e10]

                fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
                ax.hist(part_diflist,20)
                fig.savefig (self.foldername+'/diflist_hist_20_var_{0}{1}.png'.format(i,self.appendix))
                plt.close(fig)




        except: #если диод, там плоский список
            part_diflist = est['Diflist']
            fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
            ax.hist(part_diflist,20)
            fig.savefig (self.foldername+'/diflist_hist_20_{0}.png'.format(self.appendix))
            plt.close(fig)







        fig, ax = plt.subplots( nrows=1, ncols=1 )  # create figure & 1 axis
        ax.plot (est['Sklist'])
        fig.savefig (self.foldername+'/diflist_Skdrop_{0}.png'.format(self.appendix))
        plt.close(fig)

        return est

















