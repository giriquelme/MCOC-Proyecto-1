import scipy as sp
import numpy as np
import random
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

def q(t,x,r):
    if t%300==0 and r==x:
        return 0.35
    else:
        return 0.

    
    
u0 = fun_u0(x)

#Parametros del problema
dt = 1
alpha = []

#Arreglo matricial de la conductividad termica, calor especifico y la densidad de masa de diferentes materiales, especificamente del hierro, estano, acero, ladrillo refractario,Madera
#parametros=[k,c,rho]
parametros=[[79.5,450.,7800.],[64.,210.,7310.],[54.,120.,7850.],[0.8,210.,2000.],[0.13,1700.,450.]] 

#ciclo for para calcular el alpha de cada material
for i in range(len(parametros)):
    alpha.append(parametros[i][0]*dt/(parametros[i][1]*parametros[i][2]*dx**2))

contador=1
#Loop en el tiempo para cada material en especifico
for a in alpha:
    k = 0
    u_k1= u0.copy()
    u_k2= u0.copy()
    u_k3= u0.copy()
    u_k4= u0.copy()
    #--------------Condiciones de borde iniciales:-------------------- 
    #1.- Primer caso de estudio:
    #Primeros 10 datos = 5
    for i in range(10):
        u_k1[i]=5.
    #Ultimos 10 datos inician en 20, 
    for i in range(90,n+1 ):
        u_k1[i]=20.
    #2.- Segundo caso de estudio:
    u_k2[0]=5
    u_k2[n]=20
    #3.- Tercer caso de estudio:
    u_k3[0] =0
    u_k3[n] =20
    #4.- Cuarto caso de estudio:
    u_k4[0] =0 
    u_k4[n] =20
    #Creamos una copia con el objetivo de no modificar las condiciones de bordes 
    u_km1 =u_k1.copy()
    u_km2 =u_k2.copy()
    u_km3 =u_k3.copy()
    u_km4 =u_k4.copy()

    for k in range(30000):
        t = dt*k
        r=random.randint(1,100)
        #Loop en el espacio i=1 ... n  u_km1[0] = 0  u_km1[n] = 20
        for i in range(1,n):

            #Algoritmo de diferencias finitas 1D para difusion
            #Algoritmo de igualacion de pendientes con puntos vecinos
            u_km1[i] = u_k1[i] + a*(u_k1[i+1] - 2*u_k1[i] + u_k1[i-1])
            #Algoritmo de condición de borde inferior igual a 5
            u_km2[i] = u_k2[i] + a*(u_k2[i+1] - 2*u_k2[i] + u_k2[i-1])
            #Algoritmo de perdida de calor en cierto tiempo en un punto en especifico
            u_km3[i] = (u_k3[i] + a*(u_k3[i+1] - 2*u_k3[i] + u_k3[i-1]) + q(k,i,r))
            #Algoritmo de condición de borde inferior igualada a los puntos vecinos
            u_km4[i] = u_k4[i] + a*(u_k4[i+1] - 2*u_k4[i] + u_k4[i-1])
        #Primer estudio    
        #Actualizamos los puntos de los bordes igualando pendientes con puntos vecinos
        #Borde inicial
        u_km1[0] = u_km1[1] - dx*(u_km1[2]- u_km1[1])
        #Borde final
        u_km1[n] = u_km1[n-1] + dx*(u_km1[n-1] -u_km1[n-2])

        #Segundo estudio
        u_km2[n]=20

        #Tercer estudio
        u_km3[n]=20
        #Cuarto estudio
        u_km4[0] = u_km4[1] + dx*(u_km4[2]- u_km4[1])
        u_km4[n]= 20
        

        #Avanzar la solucion a k+1
        u_k1 =  u_km1
        u_k2 =  u_km2
        u_k3 =  u_km3
        u_k4 =  u_km4
        if k % 300 == 0:
            subplot (1,4,1)
            plot(x,u_k1)
            subplot (1,4,2)
            plot(x,u_k2)
           if contador==1:
                title("                  Analisis para el Hierro k = {}   t = {} \n".format(k,k*dt))
            elif contador==2:
                title("                  Analisis para el Estano k = {}   t = {} \n".format(k,k*dt))
            elif contador==3:
                title("                  Analisis para el Acero k = {}   t = {} \n".format(k,k*dt))
            elif contador==4:
                title("                  Analisis para el Ladrillo Refractario k = {}   t = {} \n".format(k,k*dt))
            else:
                titlle("                 Analisis para Madera k={}  t{}    \n".format(k,k*dt))
            subplot (1,4,3)
            plot(x,u_k3)

            subplot (1,4,4)
            plot(x,u_k4)  

    show()
    contador +=1
    
