from matplotlib.pylab import *
import numpy as np 
 
def fi(K,alpha,t):
	fii = K*(1-np.exp(-1*alpha*t))
	return fii
def Q(fi,rho,c):
	return c*rho*fi

a = 1.          #Ancho del dominio
b = 1.          #Largo del dominio
c = 1.			#Profundidad del dominio
Nx = 30         #Numero de intervalos en x
Ny = 30         #Numero de intervalos en Y
Nz = 30			#Numero de intervalor en Z
 
dx = b / Nx     #Discretizacion espacial en X
dy = a / Ny     #Discretizacion espacial en Y
dz = c / Nz		#Discretizacion espacial en Z
 
h = dx    # = dy = dz 
 
grafico1=[]
grafico2=[]
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
u_k[0,:,:] = 20.
u_k[:,0,:] = 20.
 
#Buena idea definir funciones que hagan el codigo expresivo
def printbien(u):
    print u.T[Nx::-1,:,:]
 
print u_k               #Imprime con el eje y invertido
printbien(u_k)
 
def imshowbien(u):
    imshow(u.T[Nx::-1,:,:])
    colorbar(extend='both',cmap='plasma')
    clim(10, 30)
 
#Parametros del problema (hierro)
dt = 1.0       # s
K = 79.5       # m^2 / s   
c = 450.       # J / kg C
rho = 7800.    # kg / m^3
alpha = K*dt/(c*rho*dx**2)
 
# dx =  0.166666666667
# dt = 1.0
# alpha =  0.000815384615385
 
alpha_bueno = 0.0001
dt = alpha_bueno*(c*rho*dx**2)/K
alpha = K*dt/(c*rho)
 
 
#Informar cosas interesantes
print "dt = ", dt
print "dx = ", dx
print "K = ", K
print "c = ", c
print "rho = ", rho
print "alpha = ", alpha
 
k = 0
 
# figure(1)
# imshowbien(u_k)
# title("k = {}   t = {} s".format(k, k*dt))
# savefig("movie/frame_{0:04.0f}.png".format(k))
# close(1)
 
#Loop en el tiempo 
dnext_t = 0.05   #  20.00
next_t = 0.
framenum = 0

T=5

def u_ambiente(t,T):
    return  20. + 10*np.sin((2*np.pi/T)*t)

for k in range(int32(5./dt)):
    t = dt*(k+1)
    print "k = ", k, " t = ", t
    
    u_ambiente=20+10*np.sin((2*np.pi/T)*t)
    
    #CB esencial
    u_k[0,:,:] = 0
    u_k[:,0,:] = 0
    u_k[:,:,0] = 0
    u_k[ Nx,:,:] = 0
    u_k[:, Ny,:] = 20+10*np.sin((2*np.pi/T)*t)
    u_k[:,:, Nz] = 0
    fii = fi(K,alpha,t)
    #Loop en el espacio   i = 1 ... n-1   u_km1[0] = 0  u_km1[n] = 20
    for i in range(1,Nx):
        for j in range(1,Ny):
        	for q in range(1,Nz):
	            #Algoritmo de diferencias finitas 2-D para difusion
	 
	            #Laplaciano
	            nabla_u_k = (u_k[i+1,j,q] + u_k[i-1,j,q] +u_k[i,j+1,q]+u_k[i,j-1,q] + u_k[i,j,q+1] +u_k[i,j,q-1] - 6*u_k[i,j,q])/h**2  
	 
	            #Forward euler..
	            u_km1[i,j,q] = u_k[i,j,q] + alpha*nabla_u_k + Q(fii,rho,c)
 
 	#Condiciones del caso 2


    #CB natural
    u_km1[Nx,:,:] = u_km1[Nx-1,:,:]
    u_km1[:,Ny,:] = u_km1[:,Ny-1,:]
    u_km1[:,:,Nz] = u_km1[:,:,Nz-1]
    #Avanzar la solucion a k + 1
    u_k = u_km1
 	u_k[0,:,:] = u_k[1,:,:]
    u_k[Nx,:,:] = u_k[Nx-1,:,:]
    u_k[:,0,:] = u_k[:,1,:] 
    u_k[:,Ny,:] = 20+10*np.sin((2*np.pi/T)*t)
    u_k[:,:,0] = u_k[:,:,1]
	u_k[:,:,Nz] = u_k[:,:,Nz-1]
    #agregar otros bored
 
    print "Tmax = ", u_k.max()
 
    if t > next_t:
    	grafico1.append(u_k[Nx/2,Ny/2,Nz/2])
    	grafico2.append(t)
        #figure(1)
        #imshowbien(u_k)
        #title("k = {0:4.0f}   t = {1:05.2f} s".format(k, k*dt))
        #savefig("movie{0}.png".format(framenum))
        #framenum += 1
        #next_t += dnext_t
        #close(1)
 
# figure(2)
# imshowbien(u_k)
plot(grafico1,grafico2)
# title("k = {}   t = {} s".format(k, (k+1)*dt)) 

show()

