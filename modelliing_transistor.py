__author__ = 'reiner'
#песочница для отработки функций программы, аналога planpos





from scipy import optimize

import derivations as drv


#уравнения Кирхгофа:

b=(60,60,40) #задаём вектор коэффициентов
w=(1,2) #задаём вектор значений источников питания



def make_exp_plan (wstartend1, wstartend2,  num   ):

    """делает план эксперимента, возвращает список кортежей"""
    stepw1=(wstartend1[1]-wstartend1[0])/num
    stepw2=(wstartend2[1]-wstartend2[0])/num

    for i in (wstartend1[0]+i*stepw1 for i in range(0, num)):
        print (i);

def make_exp_plan_one_generator (wstartend, num):
    stepw1=(wstartend[1]-wstartend[0])/num
    for i in range (0, num):
        yield wstartend[0]+i*stepw1

def make_exp_plan_2_generator (wstartend1, wstartend2, num):
    stepw1=(wstartend1[1]-wstartend1[0])/num
    stepw2=(wstartend2[1]-wstartend2[0])/num

    for i in range (0, num):
        yield wstartend1[0]+i*stepw1, wstartend2[0]+i*stepw2



#http://pythonworld.ru/tipy-dannyx-v-python/vse-o-funkciyax-i-ix-argumentax.html
def kirhWrapper (b):
    #эта функция возвращает функцию же.
    #Теперь надо, чтобы она возвращала функцию с подставленными аргументами. Эту функцию скормить решальнику уравнений методом Ньютона

    def kirh(y):
        return [y[0]+y[1]-y[2],
                y[0]*b[0]-y[1]*b[1]-w[0]-w[1],
                y[1]*b[1]+y[2]*b[2]+w[1]
                ]

    return kirh


def strEvaluator (funstr:list, x:list, b:list=[], c:dict={}):
    """
    Возвращает функцию, которая возвращает list результатов векторной функции, подаваемой во входных параметрах
    Сделано только для неявных функций
    funstr - вектор функций в строковом представлении
    x - вектор входных аргументов
    b - вектор коэффициентов
    c - словарь прочих постоянных
    """
    vardict = {'x':x, 'b':b}
    vardict.update(c)

    def evalfunc(yy):
        #yy - вектор выходных аргументов, которые мы находим методом Ньютона
        vdict = vardict
        vdict['y']=yy

        fff=lambda f: eval (f, None, vardict)

        return list(map(fff, funstr))

    return evalfunc #возвращаем функцию

def debilizm (r:list, let:str):
    """
    Из вектора x[...] формирует словарь {x1: .. ,...} для совместимости с drv.Jacobean
    """
    res=dict()
    for i in range(0, len(r)):
        res[let+str(i)] = r[i]
    return res


#хотели callable - получайте!
def ret_callable_jac (funstr:list, x:list, b:list=[], c:dict={}):
    def innerj (y):
        argseq=list()
        for i in range(0, len(y)):
            argseq.append('y{0}'.format(i))

        updfunstr=list(map(lambda x: x.replace('[','').replace(']',''),  funstr))
        xdict = debilizm(x,'x')
        bdict = debilizm(b,'b')
        cdict = debilizm(c,'c')
        ydict = debilizm(y,'y')
        vardict = xdict
        vardict.update(bdict)
        vardict.update(cdict)
        vardict.update(c)
        vardict.update(ydict)

        return drv.Jakobean (updfunstr, argseq, vardict)

    return innerj





def generate_uniform_plan(funcstrlist:list, xdiaptuplelist:list, b:list, c:dict, ydisps:list, nvoly=1, outfilename="", listOfOutvars=None):
    """
    Моделирует набор экспериментальных данных, получаемых по равномерному априорному плану.
    Возвращаемое значение - список словарей с ключами x, y, b, c
    Parameters:
    funcstrlist - векторная функция (список функций)
    xdiaptuplelist - вектор диапазонов, заданных кортежами


    """





    pass



def test1():
    #на этом простом примере видно, что optimize_root выдаёт верные корни при использовании в связке с strEvaluator
    #однако же при использовании его для модели транзистора резистивной есть расхождения с вариантом функции.
    #возможно, надлежит скормить якобиан однако же




    #funstr= ["x[0]+y[0]+b[0]", "x[1]+y[1]+b[1]" ]
    #x=[1,2]
    #y=[30,30]
    #b=[100,200]
    c={}

    #w преобразовано в х
    funstr= ["y[0]+y[1]-y[2]", "y[0]*b[0]-y[1]*b[1]-x[0]-x[1]", "y[1]*b[1]+y[2]*b[2]+x[1]"]
#    w=[10,60] #задаём
 #   x=w
    b=[60,60,40] #задаём вектор коэффициентов
 #   c=[]
    x=[1,2]

 #   print (ret_callable_jac(funstr, x,b,c)([10,20, 30]))




#    exit(0)


    for i in make_exp_plan_2_generator ((10,20), (60,40),  10):
        x=list(i)
        function = strEvaluator(funstr,x,b,c)
        sol = optimize.root(function, [1, 1, 1], method='lm', jac=ret_callable_jac(funstr, x,b,c))

        #TODO можно прогнать тестирование разными методами и выбрать лучший для конкретной задачи
        #sol = optimize.root(function, [1, 1, 1], method='hybr')
        print (i, sol.x, function(sol.x))


test1()


    #return strEvaluator(funstr,x,b,c)(y)





def test():
    global b


    print ("points in plan | ", " y vector | ", " left side of equation")

    for i in make_exp_plan_2_generator ((10,20), (60,40),  10):
        w=i
        #sol = optimize.root(kirh, [1, 1, 1 ], jac=jac, method=’hybr’)

        fun = kirhWrapper(b)


        sol = optimize.root(fun, [1, 1, 1 ], method='hybr')


        #можно скормить jacobean. Но он должен быть callable с подачей y на вход, очевидно


        print (i, sol.x, fun(sol.x))



        #print (list(i).extend(sol.x).extend (kirh(sol.x)))


        #написать проверялку, а точно корни-то подходящие?


#print (test1())
#test()

