__author__ = 'reiner'

import numpy as np
import scipy as scp
import random

def generlist (n):
    r=list()
    for x in range (0,n):
        r.append(random.gauss(0,1))
    return r


n=3  #размерность вектора
nvol=100 #объем выборки

#Uv=generlist(3) //покамест сгенерировали один вектор

U=list() #генерация случайного вектора U
for x in range (0,nvol):
    U.append(generlist(n))

cov_u=np.eye(n,n)  #cov matrix Ku

cov_ksi=np.array ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])

M_ksi=np.array([-10,0,10])

A=np.zeros (shape=(n,n))


#проверка корректности:


if (np.linalg.det (cov_ksi) <0 ):
    print ("Определитель матрицы меньше 0")
    exit (1)

#вычисляем А
print ()

for i in range (0, n):
    for j in range (0, i+1):
        sum=0
        for k in range (0, j-1):
            sum+=A[i,k]*A[j,k]
        if (i==j):
            A[i,j]=np.sqrt(cov_ksi[i,j]-sum)
        else:
            A[i,j]=(cov_ksi[i,j]-sum)/A[j,j]

print (A)
print ("\n\n")

ksi=list()
for i in range (0, nvol):
    ksi.append( np.dot(A, np.array(U[i]).T)+M_ksi.T)

print (M_ksi.T)
print (np.array(U[i]).T)
print (A)

print ()

#for ll in ksi:
 #   print (ll, "\n")

tseq1=list()
tseq2=list()
tseq3=list()
for ni in range (0,nvol):
    tseq1.append(ksi[ni][0])
    tseq2.append(ksi[ni][1])
    tseq3.append(ksi[ni][2])


XX=np.array([tseq1, tseq2, tseq3])
#print (XX)
COVTEST=np.cov(XX, bias=-1)

print (COVTEST)

