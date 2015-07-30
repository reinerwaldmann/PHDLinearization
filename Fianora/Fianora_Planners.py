__author__ = 'vasilev_is'


import numpy as np
import traceback
import time
import datetime
import copy

import Fianora.Fianora_static_functions as o_g #legacy


class AbstractPlanner:
    def __init__(self, N, xstart, xend, plancachefoldername='plancache', planname='plan'):
        self.plf = plancachefoldername
        self.xstart = xstart
        self.xend = xend
        self.planname = planname
        self.N = N


    def give_us_a_plan(self, nocache=False):
        """
        In normal conditions, first trying to unpack a cached plan, then trying to make plan, cache it and return
        if nocache, then just making plan without caching
        :param nocache:
        :return:
        """
        if nocache:
            return self.make_plan()

        plan = self.take_plan_from_cache()
        if not plan:
            plan = self.make_plan()
            self.write_plan_to_cache(plan)

        return plan

    def make_file_name(self):
        "makes a filename from class parameters"
        #в принципе, planname это практически всегда имя модели. В абстрактном классе, тем не менее, это отдельное поле.
        return self.plf+'/'+ self.planname+'.txt'

    def write_plan_to_cache(self, plan):
        """
         Пишет план в файл. Значения пишутся в формате python (вектор [])
        :param plan: план (список)
        :return: ничего
        """
        filename = self.make_file_name()
        if not plan:
            return
        with open (filename, 'wt') as outfile:
            for point in plan:
                outfile.write(point.__str__())
                outfile.write('\n')



    def take_plan_from_cache(self):
        """
        Читает план из файла - запись в виде столбца векторов []
        :param filename: имя файла
        :return: план (список векторов)
        """
        filename = self.make_file_name()
        with open (filename, 'rt') as infile:
            lines = infile.readlines()
            return list(map(eval, lines))


    def make_plan(self):
        raise NotImplementedError ("please implement this method")


class UniformPlanner(AbstractPlanner):

    def make_plan(self):
        """
        Создаёт равномерный план эксперимента
        :param xstart: начало диапазона x
        :param xend: конец диапазона x
        :param N: Количество точек в плане (внимание! это количество измерений каждого вектора, то есть, реальное кол-во будет N^len(xstart))
        :return: равномерный план эксперимента
        """
        res=list()

        xstartnp=np.array(self.xstart)
        xendnp=np.array(self.xend)
        xstep = list((xendnp-xstartnp)/self.N)

        evalstr="import numpy\n\n"
        lststr=""

        for i in range (len(self.xstart)):
            evalstr+="\t"*i+"for x{0} in numpy.arange(xstart[{0}], xend[{0}], xstep[{0}]):\n".format(i)
            lststr+="x{0}".format(i)+("," if i+1<len(self.xstart) else "")
        evalstr+="\t"*(i+1)+"res.append(["+lststr+"])"


        #print (evalstr)
        exec(evalstr,  locals()) #исполняет полученную сроку, собсна, формирует список входных переменных

        # for i in range(N):
        #     res.append(list(xstart+i*xstep))

        # if len(xstart)==1:
        #     for i in range (0, len(res)):
        #         res[i]=[res[i],] #костыль для диода и не нужен

        return res


class RandomPlanner(AbstractPlanner):
    def make_plan(self):
        res=list()
        for i in range(0, self.N):
            res.append(o_g.uniformVector(self.xstart, self.xend))
        return res


class DOptimalPlanner (AbstractPlanner):
    def __init__(self, N, xstart, xend, bstart, bend, model, plancachefoldername='plancache', verbose = False):
        AbstractPlanner.__init__(self, N, xstart, xend, plancachefoldername, model.name)
        self.bstart = bstart
        self.bend = bend
        self.vebose = verbose

    def __countVbForPlan(self, expplan, b):
        """
        :param expplan: план эксперимента
        :param b: b (вектор коэффициентов)
        :param b: b (вектор коэффициентов)
        :param c: словарь доп. параметров
        :param Ve: ковариационная матрица ошибок экспериментов np.array
        :param jac: функция якобиана (на входе x,b,c=None, y=None), возвращать должно np.array
        :return: значение определителя для данного плана эксперимента
        """
        G=np.zeros((len(b),len(b)))

        for point in expplan:
            jj = self.model.jacf(point,b)
            # сама модель разбирается, откуда ей выдернуть y, то есть, если он не задан,
            # None по дефолту, то она автоматически его получит

            G+=np.dot(jj.T, jj) # вставление Ve в середину вроде как особо не требуется

        try:
            return np.linalg.inv(G)
        except BaseException as e:
            print('Fatal error in countVbForPlan: ',e)
            print('b vector=',b)
            print('expplan=',expplan)
            print('G=',G)
            return None




    def __countMeanVbForAprior_S4000(self,  expplan:list):
        """


        :return: среднее значение определителя [0] и его дисперсию [1]
        """
        DS=0 #среднее определителя
        SD=0 #дисперсия определителя

        for sss in range(1, 30): #30 - количество  попыток в выборке
            b=o_g.uniformVector (self.bstart, self.bend)
            Vb=self.countVbForPlan(expplan, b)
            D=np.linalg.det(Vb)
            if D:
                DS=(D+(sss-1)*DS)/sss  #среднее определителя
                SD=((DS-D)*(DS-D)+(sss-1)*SD)/sss #дисперсия определителя

        return DS, SD

    def make_plan(self):

        """
        Реализует априорное планирование эксперимента
        :param xstart: начало диапазона x (вектор)
        :param xend: конец диапазона x  (вектор)
        :param N: размер плана (количество контрольных точек)
        :param bstart: начало диапазона b (вектор)
        :param bend: конец диапазона b (вектор)
        :param c: словарь дополнительных переменных
        :param Ve: Ковариационная матрица y, реально её диагональ (вектор)
        :param jac: Якобиан функции, принимает на вход x,b,c,y
        :param verbosePlan: если true, пишет все планы в файлы, иначе только оптимальный
        :param initplan: начальный план. По умолчанию случаен, но может быть задан (в этом случае делается одна попытка)
        :param verbose: выдавать ли в консоль информацию optimized-original
        :return: кортеж: 0: оптимизированное значение определителя Vb, 1: оптимальный план эксперимента
        """
        verbose = self.verbose
        dopt = 100000000
        planopt = None
        Ntries1 = 7

        if verbose:
            print('\n\nДанные априорного планирования:')
            print('Неоптимизированное-оптимизированное значение среднего det(Vb)')
        prevtime=time.time()

        unifplanner = UniformPlanner(self.N, self.xstart, self.xend)
        rndplanner = RandomPlanner(self.N, self.xstart, self.xend)

        plm = unifplanner

        for i in range(0,Ntries1):
            try:
                if i>0:
                    pln = rndplanner # композиция, основанная на полиморфизме - интерфейс один, но объекты подставляются разные
                    # маленькое применение паттерна "стратегия"

                m=len(self.xstart) #длина вектора входных параметров
                plan = pln.give_us_a_plan(nocache=True)

                if verbose:
                    print('plan length:', len(plan))


                unopt=self.countMeanVbForAprior_S4000(plan)[0]
                #оптимизация
                for j in range(self.N):
                    xdot=copy.deepcopy(plan[j])
                    function = lambda x: self.countMeanVbForAprior_S4000(o_g.replaceInList(plan,j,x))[0]

                    boundsarr=list()
                    for k in range(len(self.xstart)):
                        boundsarr.append((self.xstart[k],self.xend[k]))

                    #sol = optimize.minimize (function, xdot, bounds=boundsarr)
                    #plan[j]=sol.x
                    #В этом варианте не работает, хоть результат и проходит быстрее

                    plan[j]=o_g.doublesearch(self.xstart, self.xend, xdot, function)

                dcurr=self.countMeanVbForAprior_S4000(plan)[0]

                if verbose:
                    curtime = time.time() #in seconds
                    st = datetime.datetime.fromtimestamp(time.time()-prevtime).strftime('%H:%M:%S')
                    print ("{0} unoptimized-optimized: {1}   {2} time spent: {3}".format('uniform' if i==0 else '',unopt, dcurr, st))
                    prevtime=copy.copy(curtime)
                    #FIXME починить отображение времени, кривизна ж неимоверная и ниработает!!!!

                if dcurr<dopt or planopt==None:
                    dopt=dcurr
                    planopt=plan


            except BaseException as e:
                print ('This try failed, due to exception e=',e)
                tb = traceback.format_exc()
                print(tb)

        return planopt






