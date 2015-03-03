




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
import numpy as np
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

