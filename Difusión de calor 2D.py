
import scipy as sp
from matplotlib.pyplot import*
import numpy as np


#difusion 1D
#Dominio(omega)  x E(0,L=1,0)
#du/dt - k*d^2u/dt^2*cp = q(t)
# u(0,t=0)
#u(l,t)=20
#u(x,0)=u0(x)=0
# diferencias finitas(algoritmo)
#dominio(omega) --> aproximaremos  x=0 , x=dx , xn=L tendremos:
# n+1 puntos y n intervalos , con esto dx=L/n, cada uno se asocia un u de temperatura
# estructura de datos se representa por un vector de numeros reales con un arreglo de float32 con n+1
a=1.
b=1.
Nx=6
Ny=6

dx = a /Nx
dy = b / Ny

h = dx

if dx != dy:
    print "Error!!! dx != dy"
    exit(-1) # -1 le dice al SO que el programa fallo

#Funcioo de conveniencia para calcular coordenadas del punto (i,j)

#def coord(i,j):
#    return dx*i, dy*j

#i, j = 4, 2
#x, y = dx*i, dy*j
    
coord = lambda i,j : (dx*i, dy*j)
x,y =  coord(4,2)

print "x = ", x
print "y = ", y

u_k = np.zeros((Nx+1, Ny+1), dtype = np.double) # dtype es para el tipo de datos(int, float, float32, float64, double, etc)
u_km1 = np.zeros((Nx+1, Ny+1), dtype = np.double)

#CB escencial
#u_k[Ny,:] = 20
#u_k[:,0]=20

#buena idea definir funciones que hagan el codigo expresivo
def printbien(u):
    print u.T[Nx::-1,:]
    
u_k[0,:] = 20
u_k[:,0]=20

def imshowbien(u):
    imshow(u.T[Nx::-1,:])
    
print u_k
printbien(u_k)
figure()
imshowbien(u_k)
show()

#parametros del problema
dt=1    #s
K=[] #m^2/s
c=[]  #J/KgC
rho=[] #Kg/m^3
alpha=[] 
#Arreglo matricial de la conductividad termica, calor especifico y la densidad de masa de diferentes materiales, especificamente del hierro, estano, acero, ladrillo refractario, madera
#parametros=[k,c,rho]
parametros=[[79.5,450.,7800.],[64.,210.,7310.],[54.,120.,7850.],[0.8,210.,2000.],[0.13,1700.,450.]] 


#ciclo for para calcular el alpha de cada material y guardamos mediante listas k,c y rho respectivos
for i in range(len(parametros)):
    alpha.append(parametros[i][0]*dt/(parametros[i][1]*parametros[i][2]*dx**2))
    K.append(parametros[i][0])
    c.append(parametros[i][1])
    rho.append(parametros[i][2])


for l in range(len(alpha)):
    
    #Primero desprendemos los valores indices de cada lista generada
    a=alpha[l]
    print "dt=",dt
    print "dx=",dx
    print"K=",K[l]
    print "c=",c[l]
    print "rho=",rho[l]
    print "alpha=",a
    k=0
    figure(1)
    imshowbien(u_k)
    title("k= {}  t= {} s".format(k,k*dt))

    #savefig("movie/frame_ {0:04.0f}.png".format(k))
    close(1)

    for k in range (1,100):
        t=dt*(k+1)

        #loop en el espacio i=1....n   u_km1[0]=0 u_km1[n]=20
        u_k[0,:] = 20
        u_k[:,0]=20

        for i in range(1,Nx-1):
            for j in range(1, Ny-1):
                #algoritmo de dif finitas
                nabla_u_k = (u_k[i-1,j]+ u_k[i+1,j]+ u_k[i, j-1] + u_k[i,j+1] - 4*u_k[i,j])/h**2
                u_km1[i,j] = u_k[i,j] + a*nabla_u_k

        # avanzar la solucion a k+11
        u_k=u_km1

        #CB natural
        u_km1[Nx,:] = u_km1[Nx,:]   
        u_km1[:, Ny] = u_km1[:,Ny]

        figure(1)
        imshowbien(u_k)
        title("k= {}  t= {} s".format(k,k*dt))
        #Guardamos cada imagen por tipo de material, en la misma carpeta donde colocamos el archivo (.py)
        if l==0:
            savefig('Hierro{0}.png'.format(k))
        elif l==1:
            savefig('Estano{0}.png'.format(k))
        elif l==2:
            savefig('Acero{0}.png'.format(k))
        elif l==3:
            savefig('Ladrillo{0}.png'.format(k))
        elif l==4:
            savefig('Madera{0}.png'.format(k))   


    #
    #figure(2)
    #imshowbien(u_k)
    #title("k= {}  t= {} s".format(k,k*dt))
    #
    show()
    print "__________________________________"




