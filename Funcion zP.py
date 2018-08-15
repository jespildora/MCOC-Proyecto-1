#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 14 11:50:11 2018

@author: jespildora
"""
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
from numpy.linalg import inv
datos=sp.load('mck.npz')
C=datos['C']
M=datos['M']
K=datos['K']
A=datos['A']
invM=inv(M)
Cap=[150,250,500,800] # [KN]   Capacidad maxima friccional del disipador
CAP=sp.zeros(20)
CAP[:]=150
#Parametros
xi =0.025            # Razon de amortiguamiento critico
f1=0.2
f2=2.
m = 1.          # [kg]  Masa del oscilador
k = 10.         # [N/m] Rigidez del oscilador
Cap = 20.       # [N]   Capacidad maxima friccional del disipador
vr = 0.0001     # [m/s] Velocidad de referencia para la aproximacion de la friccion via tanh. 
d0 = 0         # [m]   Condicion inicial de desplazamiento
v0 = 0          # [m/s] Condicion inicial de velocidad
dt = 0.01      # [s]   Paso de integracion a usar
tmax = 100      # [s]   Tiempo maximo de integracion 
 

#Definimos la funcion del lado derecho de la EDO de primer orden
#a resolver zp = fun(t,z).  zp es la derivada temporal de z.

def fun(t,z):
    i=1
    Fr=sp.zeros(40)
    while i<20:          
          Fr[i+20]=-(CAP[i])*invM[i,i]*sp.tanh((z[20+i])-z[20+i-1])
          i+=1

    return sp.matmul(A,z)+Fr

#Vector de tiempo y su largo (Nt)
t = sp.arange(0, tmax, dt)
Nt = len(t)

#Inicializar una matriz z para guardar las solucion discretizada
z_euler = sp.zeros((40,Nt+1))
z_RK45 = sp.zeros((40,Nt+1)) 

#Condicion inciial en t = 0, i = 0. 
z0 = sp.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

 
z_euler[:,0] = z0
z_RK45[:,0] = z0
 
 
print "Integrando con Euler"
fun.tnextreport = 0
fun.solver = "Euler"
    
i = 1
ti = dt 
while (ti < tmax):
    z_euler[:,i] = dt * fun(ti, z_euler[:,i-1]) + z_euler[:,i-1]
    ti += dt
    i += 1
print 'Itegrando ',i,' iteraciones con Metodo de Euler...'
print'Integrando ',i,' iteraciones con RK45...'
fun.tnextreport = 0
fun.solver = "RK45"
solucion_rk45 = solve_ivp(fun, [0., tmax], z0, method='RK45', t_eval=t, vectorized=False )
z_RK45[:,1:] = solucion_rk45.y


#Graficar solucion en desplazamiento y velocidad para ambos metodos
plt.figure()
 
for z, lab in zip([z_euler, z_RK45], ["Euler", "RK45"]):
 
    u = z[0,:]
    v = z[20,:]
 
    #Extraer desplazamientos y velocidades
    u = z[0,:-1]
    v = z[20,:-1]
    
    umax= max(abs(u))
    plt.subplot(2,1,1)
    plt.plot(t, u, label=lab)
    plt.ylim([-1.5*umax, 1.5*umax])
    plt.xlim([0, tmax])
    plt.ylabel("Despazamiento, $u = z_1$ (m)")
    plt.grid(True)
 
    vmax = max(abs(v))
    plt.subplot(2,1,2)
    plt.plot(t, v)
    plt.ylabel("Velocidad, $\dot{u} = z_2$ (m/s)")
    plt.xlabel("Tiempo, $t$ (s)")
    plt.ylim([-1.5*vmax, 1.5*vmax])
    plt.xlim([0, tmax])
    plt.grid(True)
 
plt.subplot(2,1,1)
plt.legend()
plt.suptitle("Solucion por metodo de Euler")
 
plt.show()