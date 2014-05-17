from cmath import sqrt
import numpy as np


import matplotlib.pyplot as plt

import scipy as scp
import random

import sympy as smp


__author__ = 'reiner'


#returns u2, u3
def funct(u1, r1, r2, r3):
    return u1* (r2+r3)/(r1+r2+r3), u1* r3/(r1+r2+r3)

#print (smp.diff("u1* (r2+r3)/(r1+r2+r3)", "r2"))

def derivate(u1, r1, r2, r3):
    # rsum=(r1+r2+r3)**2
    # du3dr1 = -1*u1*r3/rsum
    # du3dr2 = -1*u1*r3/rsum
    # du3dr3 = (r1+r2)*u1/rsum
    #
    # du2dr1 = -1*u1*r2/rsum+du3dr1
    # du2dr2 = (r1+r3)*u1/rsum+du3dr2
    # du2dr3 = -1*u1*r2/rsum+du3dr3

     du2dr1 = -u1*(r2 + r3)/(r1 + r2 + r3)**2
     du2dr2 = -u1*(r2 + r3)/(r1 + r2 + r3)**2 + u1/(r1 + r2 + r3)
     du2dr3 = -u1*(r2 + r3)/(r1 + r2 + r3)**2 + u1/(r1 + r2 + r3)

     du3dr1 = -r3*u1/(r1 + r2 + r3)**2
     du3dr2 = -r3*u1/(r1 + r2 + r3)**2
     du3dr3 = -r3*u1/(r1 + r2 + r3)**2 + u1/(r1 + r2 + r3)

     return (du2dr1, du2dr2, du2dr3), (du3dr1, du3dr2, du3dr3)



#V - corr. matrix
#считает дисперсию по матожиданию и ковариационной матрице
#derivatefunc выдаёт производные см выше
def countDispLinearization (u1, r1, r2, r3, V,  derivatefunc):
    deriv=derivatefunc(u1,r1,r2,r3)
    first_memberU2=0
    for i in range (0,3):
        first_memberU2+=(deriv[0][i]**2)*(V[i,i])
    first_memberU3=0
    for i in range (0,3):
        first_memberU3+=(deriv[1][i]**2)*(V[i,i])
    sumsecMU2=0
    for i in range (0, 3):
        for j in range (0, 3):
            if (i<j):
                sumsecMU2+=deriv[0][i]*deriv[0][j]*V[i,j]
    sumsecMU3=0
    for i in range (0, 3):
        for j in range (0, 3):
            if (i<j):
                sumsecMU3+=deriv[1][i]*deriv[1][j]*V[i,j]
    return first_memberU2+sumsecMU2, first_memberU3+sumsecMU3
    #return first_memberU2, first_memberU3

#генерирует список случайных чисел м=0 д=1
def generlist (n):
    r=list()
    for x in range (0,n):
        r.append(random.gauss(0,1))
    return r


#генерирует кортеж случайных последовательностей, коррелированных по матрице и с матожиданием во вх. параметрах
def generrandvals (M_ksi, cov_ksi, nvol=1000):
    n=3  #размерность вектора
    #nvol=1000 #объем выборки
    U=list() #генерация случайного вектора U
    for x in range (0,nvol):
        U.append(generlist(n))
    cov_u=np.eye(n,n)  #cov matrix Ku

      #проверка корректности:

    if (np.linalg.det (cov_ksi) <0 ):
        print ("Определитель матрицы меньше  0")
        return None


    A=np.linalg.cholesky(cov_ksi)

    ksi=list()
    for i in range (0, nvol):
        ksi.append( np.dot(A, np.array(U[i]).T)+M_ksi.T)

    tseq1=list()
    tseq2=list()
    tseq3=list()
    for ni in range (0,nvol):
        tseq1.append(ksi[ni][0])
        tseq2.append(ksi[ni][1])
        tseq3.append(ksi[ni][2])

    XX=np.array([tseq1, tseq2, tseq3])
    COVTEST=np.cov(XX, bias=-1)

    #print ("Mean Values:")
    #print (np.mean(tseq1), np.mean(tseq2), np.mean(tseq3)  )

   # print ("Cov matrix:")
    #print (COVTEST)

    #plt.hist(tseq1, 100, label="tseq1")
    #plt.xlabel('Smarts')
    #plt.ylabel('Probability')
    #plt.title("TSEQ1")
    #plt.show()
    return tseq1, tseq2, tseq3, COVTEST, np.mean(tseq1), np.mean(tseq2), np.mean(tseq3)


#считает дисперсию и матожидание выхода, которые получаются, ежели применить сгенерированные последовательности к функции
def countDispMonteKarlo (u, rseq1, rseq2, rseq3, func):
    u2list=list()
    u3list=list()


    i=0
    for r1 in rseq1:
        r2=rseq2[i]
        r3=rseq3[i]
        i+=1
        u2list.append(func(u, r1, r2, r3)[0])
        u3list.append(func(u, r1, r2, r3)[1])

#    for r1, r2, r3 in rseq1, rseq2, rseq3:
 #       u2list.append(func(u, r1, r2, r3)[0])
 #       u3list.append(func(u, r1, r2, r3)[1])

    return np.mean(u2list), np.var(u2list), np.mean(u3list), np.var(u3list)

def countDispMonteKarlo1 (u, sectuple, func):
    return countDispMonteKarlo(u, sectuple[0], sectuple[1], sectuple[2], func)




def test1():
    print ("Monte-Karlo:")
    print ("MU2, DU2, MU3, DU3")
    rndv=generrandvals(M, V, 1000)
    print (countDispMonteKarlo1(100, rndv ,funct))
    print ("Mean Values as if lin:")
    print (funct(100, M[0], M[1], M[2]))
    print ("Linearization:")
    print ("DU2, DU3")
    print ( countDispLinearization(100, rndv[4], rndv[5], rndv[6], rndv[3], derivate))
    print ("Pretending that model is perfect:")
    print ( countDispLinearization(100, M[0], M[1], M[2], V, derivate))


#ranges:
#tseq1, tseq2, tseq3, COVTEST, np.mean(tseq1), np.mean(tseq2), np.mean(tseq3)

def test2(iM, iV, nvol):
    rndv=generrandvals(iM, iV, nvol)
    truedisp=countDispMonteKarlo1(100, rndv ,funct)
#    lineardisp=countDispLinearization(100, rndv[4], rndv[5], rndv[6], rndv[3], derivate)
    lineardisp=countDispLinearization(100, iM[0], iM[1], iM[2], iV, derivate)

    print ("true", "linear", "diff")
    print (truedisp[1], lineardisp[0], 100*(truedisp[1]-lineardisp[0])/truedisp[1], "%" )
    print (truedisp[3], lineardisp[1], 100*(truedisp[3]-lineardisp[1])/truedisp[3], "%" )


V=np.array       ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])

V1=np.array      ( [[4, 0, 0],
                    [0, 9, 0],
                    [0, 0, 16]])

#V=V1
M=np.array([20,30,400])
#test1()
#for x in [100, 1000, 10000, 100000, 1000000]:
for x in [100, 1000, 10000]:
    test2(M,V,x)





#testing sympy

#def funct(u1, r1, r2, r3):
#    return u1* (r2+r3)/(r1+r2+r3), u1* r3/(r1+r2+r3)

#print (smp.diff("u1* (r2+r3)/(r1+r2+r3)", "r2"))
#return u1* (r2+r3)/(r1+r2+r3), u1* r3/(r1+r2+r3)
#
# f1="u1* (r2+r3)/(r1+r2+r3)"
# f2="u1* r3/(r1+r2+r3)"
# r1="r1"
# r2="r2"
# r3="r3"
#
# print (smp.diff (f1, r1 )    )
# print (smp.diff (f1, r2 )    )
# print (smp.diff (f1, r3 )    )
#
# print (smp.diff (f2, r1 )    )
# print (smp.diff (f2, r2 )    )
# print (smp.diff (f2, r3 )    )
#



# for r1 in range (10,15):
#      for r2 in range (100,200,10):
#          for r3 in range (200,205):
#              M=np.array([r1, r2, r3])
#              test2(M, V, 1000)



#generrandvals(M, V)
#print (dispersion(10,20,30,40, V, derivate))

#1Теперь надо - сгенерировать последовательность r, u1 разброса не имеет
#2посчитать её ковар. матрицу (для повышения точности)
#3Все члены последовательности скормить функции func, то, что из ней вылезет, прокрутить на матожидание и дисперсию
#4Прокрутить фунцию dispersion, выдав ей матожидание и ковар. матрицу сгенерированной последовательности
#5Сравнить полученные дисперсии на шагах 3 и 4. Если они окажутся равными, значит, линеаризация в этой точке хороша.


#drawing via pyplot::
#http://matplotlib.org/1.3.1/api/pyplot_api.html




