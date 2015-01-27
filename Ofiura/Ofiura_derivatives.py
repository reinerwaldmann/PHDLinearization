__author__ = 'vasilev_is'

import sympy

def makeDerivlist (iline:str, listvars:list):
    """

    :param line:
    :param listvars: список переменных
    :return: список производных по данным переменным
    """
    line = iline.replace('[','').replace(']','').replace('math.exp','exp')

    return list (map  (lambda x: sympy.diff(line, x), listvars  ))

def makeDerivListb (line:str, b:list):
    """
    :param line:
    :param b: вектор b
    :return: список производных по данным переменным
    """
    letterlist = list()
    for i in range (len(b)):
        letterlist.append('b{0}'.format(i))

    return makeDerivlist (line, letterlist)

def makeDerivMatrix(lines:list, b:list):
    """
    :param lines: список строковых описаний функций в подготовленном виде     funlst = ["b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0","3*b0" ] #в подготовленном для взятия производной виде
    :param b:
    :return: Jacobean в символьном виде
    """
    return list (map ( lambda x: makeDerivListb(x,b),lines  ) )



#def makeDeriveMatrixPrepared (lines:list, b:list):

def backwardOpt (line):
    """
    Ищем x y b c и добавляем квадратные скобочки
    :param line:
    :return:
    """

    for i in range (100): #считаем, что уж b[100] явно не будет у нас
        pass




def test ():
    funstr = "b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0" #в подготовленном для взятия производной виде
    funlst1 = ["b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0","3*b0" ] #в подготовленном для взятия производной виде

    funlst  = ["b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0]"]


    b=[1,2,2]
    #print (makeDerivListb(funstr, b) )
    print (makeDerivMatrix(funlst, b) )
    print (makeDerivMatrix(funlst1, b) )

test()


