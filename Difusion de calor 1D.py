# -*- coding: utf-8 -*-
"""
Created on Tue Aug 06 16:28:17 2019

@author: Pedro
"""

import scipy as sp
import numpy as np
from matplotlib.pyplot import *

#Difusion 1D
#Dominio: x E {0,L(1.0)}
#Ecuacion: du/dt -k/(cp) * d^2u/dx^2 = q(t)
#Condiciones de borde:
#   u(o,t) = 0                                   Forma fuerte
#   u(L,t) = 20
#   u(x,0) = uo(x) = 0

#Diferencias finitas (algoritmos)
#Discretizacion espacial
#   dominio aproximar dx = L/n
#   x0= 0, uo
#   x1 = dx, u1 
#   xn = L, uL 

L = 1.0
n = 100
dx = L/n

#Vector con todos los x
x = np.linspace(0,L,n+1)

#Condicion inicial
def fun_u0(x):
    return 10*np.exp(-(x-0.5)**2/0.1**2)

def q(t):
    return x
    
u0 = fun_u0(x)

#creando vector de solucion u en el tiempo o paso k
u_k = u0.copy()    #Crea una nueva instancia del vector sin sobreescribirlo

#Condiciones de borde iniciales
#Primeros 10 datos = 5
for i in range(10):
    u_k[i]=5.

#Ultimos 10 datos inician en 20
for i in range(90,n):
    u_k[i]=20
u_k[n] = 20.

#Temperatura en el tiempo k +1 = dt * (k+1)
u_km1 =u_k.copy()

#Parametros del problema
dt = 1
k= 79.5
c= 450.
rho= 7800.
alpha = k*dt/(c*rho*dx**2)

plot(x,u0)

#Loop en el tiempo
k = 0
for k in range(7000):
    t = dt*k
    print "k = ", k, " t = ", t 
    
    #Loop en el espacio i=1 ... n  u_km1[0] = 0  u_km1[n] = 20
    for i in range(1,n):

        #Algoritmo de diferencias finitas 1D para difusion
        Q = q(x)
        u_km1[i] = u_k[i] + alpha*(u_k[i+1] - 2*u_k[i] + u_k[i-1])
        
    #Actualizamos los puntos de los bordes igualando pendientes con puntos vecinos
    #Borde inicial
    u_k[0] = u_km1[1] - dx*(u_km1[2]- u_km1[1])
    #Borde final
    u_k[n] = u_km1[n-1] + dx*(u_km1[n-1] -u_km1[n-2])
        
    #Avanzar la solucion a k+1
    u_k =  u_km1
    if k % 200 == 0:
        plot(x,u_k)

title("k = {}   t = {} s".format(k,k*dt))

show()