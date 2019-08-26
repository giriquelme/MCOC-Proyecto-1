from matplotlib.pylab import *
from datetime import datetime, timedelta
import numpy as np 
import csv

#Lo que hacemos aqui es definir todo lo relacionado con q(t)
def fi(K,alpha,t):
    fii = K*(1-np.exp(-1*alpha*t))
    return fii

def Q(t):
    tau = 9.8
    b = 0.98
    au = 0.7
    E = 27.1
    R = 8.31
    tr = 20
    tc = 27
    return 38591.46031*0.400463*(tau/t)**b*(b/t)*au*np.exp(-(tau/t)**b)*np.exp(E/R*(1/(273+tr)- 1/(273+tc)))
    
# Funcion que nos entrega el grado de hidratacion con respecto al tiempo con un grado de hidratacion maximo de 0.7
def H(t):
    tau = 9.8
    b = 0.98
    au = 0.7    
    return au*np.exp(-(tau/t)**b)

# Funcion que nos entrega el calor especifico del hormigon en funcion del grado de hidratacion y de la temperatura
def C(h, T):   # h es el grado de hidratacion del hormigon y T es la temperatura
    Wc = 498.06
    Wa = 738.23
    Ww = 7.42
    A = 8.4
    B = 339
    Cc = 1140
    Ca = 770
    Cw = 4186
    densidad = 2400
    return (Wc*h*(A*T+B) + Wc*(1-h)*Cc + Wa*Ca + Ww*Cw)/densidad
    
#Buena idea definir funciones que hagan el codigo expresivo
def printbien(u):
    print u.T[Nx::-1,:,:]
 
def imshowbien(u):
    imshow(u.T[Nx::-1,:,:])
    colorbar(extend='both',cmap='plasma')
    clim(10, 30)

TAmbiente=[]
with open('TemperaturaAmbiente.csv') as File:
    reader = csv.reader(File)
    tiempoInicial= datetime(2018,11,27,10,13,8)
    contador = 0
    for row in reader:
        if contador != 0:
            separator1 = row[0].split("/")
            separator2 = row[1].split(":")
            fecha = datetime( 2018 ,int32(separator1[1]),int32(separator1[0]),int32(separator2[0]),int32(separator2[1]),int32(separator2[2]))
            convercion =  fecha - tiempoInicial
            segundos = convercion.total_seconds()
            TAmbiente.append([float(row[3]),segundos])
        contador +=1
print TAmbiente

sTotales=1136375
a = 1.          #Ancho del cubo de Hormigon
b = 0.5          #Largo del cubo de Hormigon
c = 0.5         #Profundidad del cubo de Hormigon (redondeado por temas de no complicar mucho el codigo)
Nx = 50       #Numero de intervalos en x
Ny = 25        #Numero de intervalos en Y
Nz = 25     #Numero de intervalor en Z
 
dx = a / Nx     #Discretizacion espacial en X
dy = b / Ny     #Discretizacion espacial en Y
dz = c / Nz     #Discretizacion espacial en Z
 
h = dx    # = dy = dz 
 
if dx != dy or dx != dz or dz != dy :
    print("ERRROR!!!!! dx != dy  or dx != dz or dz != dy")
    exit(-1)   #-1 le dice al SO que el programa fallo.....

#Funcion de conveniencia para calcular coordenadas del punto (i,j)

# def coords(i,j):
#   return dx*i, dy*j
# x, y = coords(4,2)  
 
# i, j = 4, 2 
# x, y = dx*i, dy*j

coords = lambda i, j, q : (dx*i, dy*j, dz*q)
x, y, z = coords(4,2,6) 

print "x = ", x
print "y = ", y
print "z = ", z
 
u_k = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
u_km1 = zeros((Nx+1,Ny+1,Nz+1), dtype=double)   #dtype es el tipo de datos (double, float, int32, int16...)
 
#CB esencial
u_k[:,:,:] = 27.
#u_k[0,:,:] = 20.
#u_k[-1,:,:] = 20.
#u_k[:,0,:] = 20.
#u_k[:,-1,:] = 20.
#u_k[:,:,0] = 20.

#Condicion inicial

#Parametros del problema (hierro)
dt = 1.0       # s
K = 78.8     # m^2 / s   
#c = 450.       # J / kg C
rho = 2400.    # kg / m^3
#alpha = K*dt/(c*rho)

# dx =  0.166666666667
# dt = 1.0
# alpha =  0.000815384615385
 
alpha_bueno = 0.0001
#dt = alpha_bueno*(c*rho*dx**2)/K
#alpha = K*dt/(c*rho)
 
 
#Informar cosas interesantes
print "dt = ", dt
print "dx = ", dx
print "K = ", K
#print "c = ", c
print "rho = ", rho
#print "alpha = ", alpha
 
k = 0
 
# figure(1)
# imshowbien(u_k)
# title("k = {}   t = {} s".format(k, k*dt))
# savefig("movie/frame_{0:04.0f}.png".format(k))
# close(1)
 
#Loop en el tiempo 
dnext_t = 20.  #  20.00
next_t = 0.
framenum = 0

puntomedio1=[]
puntomedio2=[]
puntomedio3=[]
puntocero1=[]
puntocero2=[]
puntocero3=[]
punto1=[]       
punto2=[]
punto3=[]
tiempo=[]



for k in range(len(TAmbiente)):
    t = (TAmbiente[k][1] +1)/60
    dt = (TAmbiente[k][1] - TAmbiente[k][1])/60
    print "k = ", k, " t = ", t
    u_k[:,:, Nz] = TAmbiente[k][0]
    #Loop en el espacio   i = 1 ... n-1   u_km1[0] = 0  u_km1[n] = 20
    for q in range(1,Nz):
        for j in range(1,Ny):
            for i in range(1,Nx):
                #Calculamos alpha con el calor especifico respectivo al grado de hidratacion del hormigon
                alpha = K*dt/(C(H(t), u_k[i,j,q])*rho*dx**2)   
                #Algoritmo de diferencias finitas 3-D para difusion
                #Laplaciano
                nabla_u_k = (u_k[i+1,j,q] + u_k[i-1,j,q] +u_k[i,j+1,q]+u_k[i,j-1,q] + u_k[i,j,q+1] +u_k[i,j,q-1] - 6*u_k[i,j,q])/h**2  
                #Forward euler..
                u_km1[i,j,q] = u_k[i,j,q] + alpha*nabla_u_k + Q(t)
                
    #Condiciones del caso 2
    #u_k[0,:,:] = 20.
    #u_k[-1,:,:] = 20.
    #u_k[:,0,:] = 20.
    #u_k[:,-1,:] = 20.
    #u_k[:,:,0] = 20.

    #CB natural
    u_km1[0,:,:] = u_km1[1,:,:]
    u_km1[-1,:,:] = u_km1[-1,:,:]
    u_km1[:,0,:] = u_km1[:,1,:]
    u_km1[:,-1,:] = u_km1[:,-2,:]
    u_km1[:,:,0] = u_km1[:,:,1]
    
    #Avanzar la solucion a k + 1
    u_k = u_km1
    #agregar otros bored
 
    print "Tmax = ", u_k.max()
    
    puntomedio1.append(u_k[Nx/2,Ny/2,Nz/2])
    puntomedio2.append(u_k[Nx/2,Ny/2,2])
    puntomedio3.append(u_k[Nx/2,Ny/2,Nz-2])

    puntocero1.append(u_k[2,2,2])
    puntocero2.append(u_k[2,2,Nz/2])
    puntocero3.append(u_k[2,2,Nz-2])

    punto1.append(u_k[Nx/2,2,2])
    punto2.append(u_k[Nx/2,2,Nz/2])
    punto3.append(u_k[Nx/2,2,Nz-2])
    tiempo.append(t)
        #figure(1)
        #imshowbien(u_k)
        #title("k = {0:4.0f}   t = {1:05.2f} s".format(k, k*dt))
        #savefig("movie{0}.png".format(framenum))
        #framenum += 1
    next_t += dnext_t
        #close(1) 
 
#    if t > next_t:

    
    if TAmbiente[k][1] > 3000:
        break
    
# figure(2)
# imshowbien(u_k)
plot(tiempo,puntomedio1,'b')
plot(tiempo,puntomedio2,'g')
plot(tiempo,puntomedio3,'r')
plot(tiempo,puntocero1,'c')
plot(tiempo,puntocero2,'m')
plot(tiempo,puntocero3,'y')
plot(tiempo,punto1,'k')
plot(tiempo,punto2,'gray')
plot(tiempo,punto3,'violet')

title("HORMIGONES MASIVOS k = {}   t = {} s".format(k, TAmbiente[3000][1])) 

show()
