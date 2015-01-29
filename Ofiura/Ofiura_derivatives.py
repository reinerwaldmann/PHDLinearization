__author__ = 'vasilev_is'

import sympy

def makeDerivlist (iline:str, listvars:list):
    """

    :param line:
    :param listvars: список переменных
    :return: список производных по данным переменным
    """
    line = iline.replace('[','').replace(']','').replace('math.','')

    # print (line)
    # print (listvars)

    res= list (map  (lambda x: str(sympy.diff(line, x)), listvars  ))



    return list(map(backwardOpt, res))

def makeDerivListb (line:str, b:list, letter='b'):
    """
    :param line:
    :param b: вектор b
    :return: список производных по данным переменным
    """
    letterlist = list()
    for i in range (len(b)):
        letterlist.append('{0}{1}'.format(letter, i))

    return makeDerivlist (line, letterlist)

def makeDerivMatrix(lines:list, b:list, letter='b'):
    """
    :param lines: список строковых описаний функций в подготовленном виде     funlst = ["b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0","3*b0" ] #в подготовленном для взятия производной виде
    :param b:
    :return: Jacobean в символьном виде
    """
    res = list (map ( lambda x: makeDerivListb(x,b, letter),lines  ) )
    return res.__str__().replace('\'','')



#def makeDeriveMatrixPrepared (lines:list, b:list):

def backwardOpt (line):
    """
    Ищем x y b c и добавляем квадратные скобочки
    :param line:
    :return:
    """

    for i in range (100): #считаем, что уж b[100] явно не будет у нас
        line=line.replace ('x{0}'.format(i),'x[{0}]'.format(i))
        line=line.replace ('b{0}'.format(i),'b[{0}]'.format(i))
        line=line.replace ('c{0}'.format(i),'c[{0}]'.format(i))
        line=line.replace ('y{0}'.format(i),'y[{0}]'.format(i))
    return line.replace('exp', 'math.exp')





def test ():
    #funstr = "b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0" #в подготовленном для взятия производной виде
    #funlst1 = ["b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0","3*b0" ] #в подготовленном для взятия производной виде
    #funlst  = ["b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0]"]
    funlst=["b[0]*(math.exp((x[0]-y[0]*b[2])/(FT*b[1])) -1)-y[0]"]
    b=[1,2,2]
    y=[1]
    #print (makeDerivListb(funstr,
    # b) )
    print (makeDerivMatrix(funlst, y, 'y'))
    funstr = "b0*(exp((x0-y0*b2)/(FT*b1)) -1)-y0"
    print (sympy.diff(funstr,'y0').__str__())







