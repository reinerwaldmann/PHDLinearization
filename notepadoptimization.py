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

median=10 #истинное значение
relError=0.1 #доля ошибки, то есть 10%






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