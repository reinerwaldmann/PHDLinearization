

__author__ = 'vasilev_is'

# Пусть тебе надо ряд из t величин с матожиданием (m1,m2,...,mt) и ковариационной
# матрицей V[t x t] (она должна быть неотрицательно определенной и симметричной).
#
# Генеришь ряд X=(X1,X2,...,Xt) из независимых случайных с нулевым матожиданием
# и единичной дисперсией. Находишь матрицу W[t x t] такую, что W'W=V (называется это разложением Холецкого), и берешь Y = XW + (m1,m2,...,mt).

import numpy as np
import random
import numpy.random as nprnd

#сделали ряд
def generlist (n):
    r=list()
    for x in range (0,n):
        r.append(random.gauss(0,1))
    return r

# V=np.array([[ 1.09031615, -0.2378683,   0.54322079],
#  [-0.2378683,   0.29493523, -0.26074706],
#  [ 0.54322079, -0.26074706,  1.1808071 ]])


V=np.array ([[4, 2, 3],
             [2, 9, 6],
             3, 6, 16])

print ("Тестовые массивы")

num=20000

seq1 = generlist(num)
seq2 = generlist(num)
seq3 = generlist(num)
X=np.array([seq1, seq2, seq3])
#print (X)
COV=np.cov(X , bias=-1)
print (COV)
print ("\n\n\n")


#print (np.linalg.det(V)) #проверяем, что матрица неотрицательно определена
W = np.linalg.cholesky(COV) #генерируем разложение Холецкого
print ("Cholesky")
print (W)
# и единичной дисперсией. Находишь матрицу W[t x t] такую, что W'W=V (называется это разложением Холецкого), и берешь Y = XW + (m1,m2,...,mt).

tseq1=list()
tseq2=list()
tseq3=list()
for x in range (0,num):
    X=np.array(generlist(3))
    Y= np.dot(X ,  W)  #АХТУНГ!!!! X*W - это поэлементно перемножение!!!! ОЛОЛ1111!!!! ПЫЩЬ!!!!!
    tseq1.append(Y[0])
    tseq2.append(Y[1])
    tseq3.append(Y[2])


X=np.array([tseq1, tseq2, tseq3])
#print (X)
COVTEST=np.cov(X, bias=-1)
print (COVTEST)

print ("\n\n")

print (COV-COVTEST)



 #   ____  __  __  _____ _  __
 #  / __ \|  \/  |/ ____| |/ /
 # | |  | | \  / | (___ | ' /
 # | |  | | |\/| |\___ \|  <
 # | |__| | |  | |____) | . \
 #  \____/|_|  |_|_____/|_|\_\






