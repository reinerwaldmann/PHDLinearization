X=2
X = X * X / 9
J = (((((0.00021 * X - 0.0039444) * X + 0.0444479) * X - 0.3163866) * X + 1.2656208) * X - 2.2499997) * X + 1
print (J)
print (jv(0,2))



XX=1
X=XX
X = X * X / 9
J = XX * ((((((0.00001109 * X - 0.00031761) * X + 0.00443319) * X - 0.03954289) * X + 0.21093573) * X - 0.56249985) * X + 0.5)
print (J)
print (jv(1,1))








exit(0)
print ([1,2,3,4,5,6][-1])

exit(0)

mpm.mp.dps = 5
#mpm.mp.prec = 50

#mpm.pretty = True

a = mpm.mpf('1.123456789123456789123456789123456789')

b = mpm.mpf("1.123456789123456789123456789123456789")


print (a)

exit(0)
import seaborn as sns


def main():
    #for num in [10, 50, 100, 1000]:
    for num in [10]:
        # Generate some data
        x = np.random.normal(0, 0.5, num-3)

        # Add three outliers...
        x = np.r_[x, -3, -10, 12]
        plot(x)

    plt.show()

def mad_based_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.

        https://en.wikipedia.org/wiki/Standard_score#/media/File:Normal_distribution_and_scales.gif

        http://blog.caseystella.com/pyspark-openpayments-analysis-part-4.html


    """

    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)

    print ((points - median)**2)
    print (np.sum((points - median)**2 , axis=-1))

    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)

    print (diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh

def percentile_based_outlier(data, threshold=95):
    diff = (100 - threshold) / 2.0
    minval, maxval = np.percentile(data, [diff, 100 - diff])



    return (data < minval) | (data > maxval)

def plot(x):
    fig, axes = plt.subplots(nrows=2)
    for ax, func in zip(axes, [percentile_based_outlier, mad_based_outlier]):
        sns.distplot(x, ax=ax, rug=True, hist=False)

        outliers = x[func(x)] # func(x) вернуло маску, то есть True True False, притом  True на тех элементах
        # которые выбиваются. Потом этим оборотом x[func(x)] мы извлекли те элементы, которые тру.

        ax.plot(outliers, np.zeros_like(outliers), 'ro', clip_on=False)

    kwargs = dict(y=0.95, x=0.05, ha='left', va='top')
    axes[0].set_title('Percentile-based Outliers', **kwargs)
    axes[1].set_title('MAD-based Outliers', **kwargs)
    fig.suptitle('Comparing Outlier Tests with n={}'.format(len(x)), size=14)

main()


exit(0)

a=.1
b=1000

import time
import os
#ОСТЕОХОНДРОЗНОЕ


def clear():
    os.system('clear')

for k in range (10):

    for i in range (10):
        time.sleep(1)
        i += 1
        os.system('beep -f 1000')

    for i in range (10):
        time.sleep(1)
        i += 1
        os.system('beep -f 200')


    print ('Numpers of approaches:',k)
print ('FINISHED')








exit(0)

filename = 'temp1.html'

with open(filename, 'w') as f:
    pass
#   зачистка



sys.stdout = open(filename, 'a') # Перенаправить вывод в файл


cells_per_row = 5

print ('<html>')

print ("<table style='page-break-after: always'>")
print ('<tr>')
for i in range(101):

    if not i%cells_per_row:
        print ('</tr>')
        print ('<tr>')

    if not i%(23*cells_per_row):
        print ("<table style='page-break-after: always'>")
        print ("<table>")

    print ("<td style='border-style:solid; border-width:1px; padding:10px; font-weight: bold'>")
    #print ('23-{0}'.format(str(i).zfill(5)))
    print ('23-ОКР-{0}'.format(str(i)))



    print ('</td>')




print ('</tr>')
print ('</table>')
print ('</html>')














exit(0)
print(((780/2)**2+100**2)**.5)





exit(0)

#великолепно же!!! (стянуто с лурка)


import sys;

class Cout(object):
    def __lshift__(self, anything):
        sys.stdout.write(str(anything));
        return self;

cout = Cout();
endl = '\n';

cout << 'Hello, World!' << endl;

gg=[2,3,4]
cout <<gg



exit(0)












import mpmath as mpm





mpm.mp.dps=50

gf=mpm.mpf('1.0000000000000000  001')
print (gf)



a =mpm.mpf (9.999999999)
b= mpm.mpf (0.0000000001)

print (a, b)

c=a+b
print (c)




deltam=1
with mpm.workdps(mpm.mp.dps+deltam):
    c=a+b
    print (c)

print (c)










#
# with mpm.workdps(30):
#     print (g)
#
# with mpm.extradps(10):
#     print (g)
#
# print (g)








exit(0)


import matplotlib.pyplot as plt

#x = range(100)
#y = [sqrt(i) for i in x]

x=[0,5]
y=[10/3,0]

s=0
d=None
for xx in np.arange (0,5,.001):
   yy=(1-xx/5)*3.3
   if xx*yy>s:
       s=xx*yy
       d=(xx,yy)
print (d)




plt.plot(x,y,color='k')

plt.plot ([0, d[0], d[0]], [d[1], d[1], 0 ], color='r' )


plt.axis([0,5,0,5])
plt.fill_between(x,y,0,color='0.8')
plt.show()



# import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
# import numpy as np
# import matplotlib.pyplot as plt
#
# mpl.rcParams['legend.fontsize'] = 10
#
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# theta = np.linspace(-4 * np.pi, 4 * np.pi, 100)
# z = np.linspace(-2, 2, 100)
# r = z**2 + 1
# x = r * np.sin(theta)
# y = r * np.cos(theta)
# ax.plot(x, y, z, label='parametric curve')
# ax.legend()
#
# plt.show()



exit(0)


g = {'one':1, 'two':2, 'three':3, '1':23}

lst = g.keys()
[print (i) for i in lst if 't' in i] #печатает любой iterable красиво - каждый элемент на своей строчке
#красивенная конструкция, которая фильтрует и печатает
exit(0)




klist = sorted(list(g.keys()))

rslist = [g[i] for i in sorted(g.keys())]


# rslist=list()
# for i in klist:
#     rslist.append(g[i])

print (klist)

print (rslist)

num=2
print ([klist[i::num] for i in range(num)])






exit(0)


from numpy import arange, sin, pi, random


x = arange(0, 6e-2, 6e-2 / 30)
A, k, theta =  10, 1.0 / 3e-2, pi / 2
y_true = A * sin(2 * pi * k * x + theta)
y_meas = y_true + 2*random.randn(len(x))



def residuals(p, y, x):
    A, k, theta = p
    err = y - A * sin(2 * pi * k * x + theta)
    return err



def peval(x, p):
    return p[0] * sin(2 * pi * p[1] * x + p[2])


# p0 = [8, 1 / 2.3e-2, pi / 3]
# #print(array(p0))
#
#
# plsq = leastsq(residuals, p0, args=(y_meas, x))
# print(plsq[0])
#
# print(array([A, k, theta]))
#
#моделируем погрешность измерения

# median=10 #истинное значение
# relError=0.1 #доля ошибки, то есть 10%
#


sum1=0  #первая формула
sum2=0 #вторая формула

# pa=0.9999
# pb=0.9997

pa=0.8
pb=0.4

# c=0
# N=0

c=7
N=17

#
# def c_in (i,n):
#     return math.factorial(n) / (math.factorial(i) * math.factorial(n-i))
#
#
# for i in range (c-1):
#
#     sum1+=c_in(i,N)*math.pow(pa,N-i)*(math.pow(1-pa,i))
#     sum2+=c_in(i,N)*math.pow(pb,N-i)*(math.pow(1-pb,i))
#
#
#
# print (sum1, sum2)
# print (1-0.2, 0.2)
#
import decimal as dc



np.float64

mm = np.array ([1.0,2.0,np.longdouble(0.00000003+1)], dtype=np.dtype(np.longdouble))
print (mm)
print (type(mm[2]))


mm = np.array ([1.0,2.0,0.00000003+1])
print (mm)
print (type(mm[2]))


ff=dc.Decimal (0.00000000000003)
print (ff+1)

ma=np.array([ff])


print (ma+ma-ma-ma+np.array([1]))
print (type(ma[0]))




import numpy as np
im

sbarr=np.array ([1,0.00000001])
sbarr1=np.array ([1,1])

print (type(sbarr1[1]))

sbarr+=sbarr1

print (sbarr)

print (type(sbarr[1]))






exit(0)


#!/usr/bin/env python


TCP_IP = '127.0.0.1'
TCP_PORT = 5005
BUFFER_SIZE = 1024
MESSAGE = "Hello, World!"

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((TCP_IP, TCP_PORT))

s.send(bytes(MESSAGE, 'UTF-8'))





data = s.recv(BUFFER_SIZE)
s.close()

print ("received data:", data)


#Server part

__author__ = 'vasilev_is'

import socket


TCP_IP = '127.0.0.1'
TCP_PORT = 5005
BUFFER_SIZE = 20  # Normally 1024, but we want fast response

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind((TCP_IP, TCP_PORT))
s.listen(1)

conn, addr = s.accept()
print ('Connection address:', addr)
while 1:
    data = conn.recv(BUFFER_SIZE)
    if not data: break
    print ("received data:", data)
    conn.send(data[1:5])  # echo
conn.close()

