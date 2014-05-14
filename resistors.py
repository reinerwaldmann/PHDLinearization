from cmath import sqrt
import numpy as np
import scipy as scp
import random


__author__ = 'reiner'


#returns u2, u3
def funct(u1, r1, r2, r3):
    return u1* (r2+r3)/(r1+r2+r3), u1* r3/(r1+r2+r3)


def derivate(u1, r1, r2, r3):
    rsum=(r1+r2+r3)**2
    du3dr1 = -1*u1*r3/rsum
    du3dr2 = -1*u1*r3/rsum
    du3dr3 = (r1+r2)*u1/rsum

    du2dr1 = -1*u1*r2/rsum+du3dr1
    du2dr2 = (r1+r3)*u1/rsum+du3dr2
    du2dr3 = -1*u1*r2/rsum+du3dr3
    return (du2dr1, du2dr2, du2dr3), (du3dr1, du3dr2, du3dr3)



#V - corr. matrix
def dispersion (u1, r1, r2, r3, V,  derivatefunc, funct):
    deriv=derivatefunc(u1,r1,r2,r3)

    lst= list (map(lambda x: x**2, deriv[0]))
    dispssqrt=[V[i,i]**2 for i in range (0, V.ndim+1)]





  #sssss
    first_member = sqrt ( sum (np.array(deriv[0])* np.array (dispssqrt)))

    return first_member


V=np.array ( [[4, 2, 3],
                    [2, 9, 6],
                    [3, 6, 16]])


print (dispersion(10,20,30,40, V, derivate, funct))



