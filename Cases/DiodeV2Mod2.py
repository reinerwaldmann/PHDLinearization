__author__ = 'vasilev_is'

import math

import numpy as np
from scipy import optimize
import sympy

#Part1: прямая ветвь ВАХ диода
#Стандарт: implicit
#Таблица переменных:



#Описание стандартных функций приводятся в Методике программирования оценочных скриптов
#Стандартные глобальные переменные:
numnone=0 #количество раз, когда функция вернула None, не справившись с оценкой тока


#Нестандартные глобальные переменные:
FT=0.02586419 #подогнанное по pspice




def func_Kirch_DiodeV2Mod2DirectBranch(y,x,b,c):
    """
    [Реестровая]
    Уравнение Кирхгофа
    :param y:
    :param x:
    :param b:
    :param c:
    :return:
    """
    global FT
    #mm=float(b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0])
    In =   float(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )-y[0]) #+
    King = float(math.sqrt (b[2]/(b[2]+In))) #+
    Irec = float(b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )-y[0]) #+
    Kgen = float(((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2) ) #+
    Ifwd = float(In*King+Irec*Kgen) #+
    Ifwd_direct=(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )-y[0]) * math.sqrt (b[2]/(b[2]+(b[0]*(math.exp( (x[0]-y[0]*b[7])/(FT*b[1]))-1 )-y[0])))+\
                (b[3]*(math.exp((x[0]-y[0]*b[7])/(FT*b[4]))-1 )-y[0]) * (((1- (x[0]-y[0]*b[7])/b[5])**2+0.005 )**(b[6]/2))
    return [Ifwd]




def solver_Kirch_DiodeV2Mod2DirectBranch (x,b,c=None):
    """
    [Реестровая]
    :param x:
    :param b:
    :param c:
    :return:
    """
    global numnone

    global FT

    dfdy=lambda y,x,b,c=None: np.array ([[ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]])

    func = func_Kirch_DiodeV2Mod2DirectBranch
    solvinit=[1]

    try:
        solx=optimize.root(func, solvinit, args=(x,b,c), jac=dfdy, method='lm').x
        #http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
        #http://codereview.stackexchange.com/questions/28207/is-this-the-fastest-way-to-find-the-closest-point-to-a-list-of-points-using-nump
    except BaseException as e:
        #print ('diodeResistorIMPLICITfunction: Error in findroot=',e)
        numnone+=1
        return None

    if solx-solvinit==[0]*len(solx):
         numnone+=1
         return None

    return solx



def Jac_Kirch_DiodeV2Mod2DirectBranch (x,b,c,y):
    global FT


    dfdb=lambda y,x,b,c: np.matrix( [ [math.exp((-b[2]*y[0] + x[0])/(FT*b[1])) - 1,
                                      -b[0]*(-b[2]*y[0] + x[0])*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1]**2),
                                      -b[0]*y[0]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])]     ])


    dfdy=lambda y,x,b,c: np.matrix ([ -1 - b[0]*b[2]*math.exp((-b[2]*y[0] + x[0])/(FT*b[1]))/(FT*b[1])])

    #возвращает структурную матрицу
    #jacf=lambda x,b,c,y: jjacf(x,b,c,y,dfdb,dfdy)


    jacf=np.dot(np.linalg.inv(dfdy(y,x,b,c)), dfdb(y,x,b,c) )


    return jacf


def test_Kirch_DiodeV2Mod2DirectBranch():

    b=[1.238e-14, 1.8, 1000]
    #0.0026638081177255196
    rng=np.arange(0.01,1.5,0.01)
    #снимем ВАХ
    resrng=[diodeResistorIMPLICITfunction ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
#    resrngorig=[casesDiode.diode([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
    b=[1.238e-14, 1.8, 3000]
    resrng1=[diodeResistorIMPLICITfunction ([x],b)[0] for x in rng] # изменяем напряжение на базе при постоянном напряжении на коллекторе - снимаем ток базы.
 #   plt.plot(rng , resrngorig, label='r=0')
    plt.plot(rng , resrng, label='r=1000')
    plt.plot(rng , resrng1, label='r=3000')
    plt.legend(loc='upper left')
    #plt.axis([0.0,1.0,0,5])
    plt.grid()
    plt.show()

def extraction_Kirch_DiodeV2Mod2DirectBranch():
    pass



funstr="x**4"

print (str(sympy.diff(funstr, 'x')))
