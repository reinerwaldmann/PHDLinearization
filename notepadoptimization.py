#http://docs.scipy.org/doc/scipy/reference/tutorial/optimize.html

from numpy import arange, sin, pi, random, array

from scipy.optimize import leastsq

import ApriorPlanning as ap


x = arange(0, 6e-2, 6e-2 / 30)
A, k, theta = 10, 0.02,  6
y_true = A * sin(2 * pi * k * x + theta)
y_meas = y_true + 2*random.randn(len(x))

def residuals(p, y, x):
     A, k, theta = p
     err = y - A * sin(2 * pi * k * x + theta)
     return err

def peval(x, p):
     return p[0] * sin(2 * pi * p[1] * x + p[2])



def funcf (x,b):
    return array ([x[0]*(b[0]+b[1]), x[0]*b[1]])

def residuals_mine(b,x,y):
    """
    :param b: параметры, которые мы подбираем
    :param x: значения величин (список векторов)
    :param y: реальные данные (список векторов)
    :return:
    """



#http://muzhig.ru/least-squares-fitting-bell-curves/


p0 = [8, 1 / 2.3e-2, pi / 3]
#print(array(p0))
plsq = leastsq(residuals, p0, args=(y_meas, x))
#print(plsq[0])

#print(array([A, k, theta]))

print (x)
print (y_meas)









jacf = lambda x,b, y: np.matrix([ [x[0], x[0]], [0, x[0] ]   ])



print  (grandCountGN_Ultra(funcf, jacf, expdata, [1,1]))











# import matplotlib.pyplot as plt
# plt.plot(x, peval(x, plsq[0]),x,y_meas,'o',x,y_true)
# plt.title('Least-squares fit to noisy data')
# plt.legend(['Fit', 'Noisy', 'True'])
# plt.show()
#
