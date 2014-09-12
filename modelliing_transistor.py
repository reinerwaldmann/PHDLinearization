__author__ = 'reiner'
#программая является аналогом kirhgres


from scipy import optimize
import derivations as deriv

#уравнения Кирхгофа:

b=(60,60,40) #задаём вектор коэффициентов
w=(1,2) #задаём вектор значений источников питания

def kirh(y):
    global b

    return [y[0]+y[1]-y[2],
            y[0]*b[0]-y[1]*b[1]-w[0]-w[1],
            y[1]*b[1]+y[2]*b[2]+w[1]
            ]

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



def test():
    print ("points in plan | ", " y vector | ", " left side of equation")

    for i in make_exp_plan_2_generator ((10,20), (60,40),  10):
        w=i
        #sol = optimize.root(kirh, [1, 1, 1 ], jac=jac, method=’hybr’)
        sol = optimize.root(kirh, [1, 1, 1 ], method='hybr')

        print (i, sol.x, kirh(sol.x))


        #print (list(i).extend(sol.x).extend (kirh(sol.x)))


        #написать проверялку, а точно корни-то подходящие?



test()

