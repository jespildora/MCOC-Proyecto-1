#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 12 15:21:37 2018

@author: jespildora
"""
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import numpy as np
from numpy.linalg import inv
import math
e=23.5 #Modulo de Young en GPa

#Inercias
i1=1080000.0 #cm4 60x60
i2=2000833.3 #cm4 70x70
i3=3413333.3 #cm4 80x80
i4=5467500.0 #cm4 90x90
i5=8333333.3 #cm4 100x100

ei1=e*i1/100000 #KN*m2
ei2=e*i2/100000 #KN*m2
ei3=e*i3/100000 #KN*m2
ei4=e*i4/100000 #KN*m2
ei5=e*i5/100000 #KN*m2

l1=4.0 #m Longitud de la columna
l2=2.8 #m LOngitud de la columna

Cap=[150,250,500,800] #Expresadas en KN
#Luego se modificaran las unidades de las cpaacidades para que tengan las mismas que K
CapMax=[5000] #KN

#Calculamos 10 rigideces distintas mediante (12EI/l**3), las primeras 5 corresponden a las columnas de largo 4m
#con sus respectivas areas transversales, desde 60x60 hasta 100x100
#Las siguientes 5, sigue el mismo patron pero con columnas de largo 2,8m
listak=[12*ei1/(l1**3),12*ei2/(l1**3),12*ei3/(l1**3),12*ei4/(l1**3),12*ei5/(l1**3),12*ei1/(l2**3),12*ei2/(l2**3),12*ei3/(l2**3),12*ei4/(l2**3),12*ei5/(l2**3)]


#Se calculan los valores de K, que es una matriz tri-diagonal, se calcularon por piso
#Sumando las rigidecez que afectaria el movimiento de los pisos

#Valores de K en MN/m , 'Mega newton por metro'
K = np.zeros((20,20))
K[0,0]=6*listak[0]+2*listak[1]+4*listak[4]+listak[3]+6*listak[5]+2*listak[6]+4*listak[9]+listak[8]
K[0,1]=-(6*listak[5]+2*listak[6]+4*listak[9]+listak[8])
K[1,0]=-(6*listak[5]+2*listak[6]+4*listak[9]+listak[8])
i=1
while i<=2:
    K[i,i]=12*listak[5]+4*listak[6]+8*listak[9]+2*listak[8]
    K[i,i+1]=-(6*listak[5]+2*listak[6]+4*listak[9]+listak[8])
    K[i+1,i]=-(6*listak[5]+2*listak[6]+4*listak[9]+listak[8])
    i+=1
K[3,3]=10*listak[5]+2*listak[6]+8*listak[9]+listak[7]+listak[8]
K[3,4]=-(4*listak[5]+4*listak[9]+listak[7])
K[4,3]=-(4*listak[5]+4*listak[9]+listak[7])
i+=1
while i<=6:
    K[i,i]=8*listak[5]+8*listak[9]+2*listak[7]
    i+=1
K[4,5]=-(4*listak[5]+4*listak[9]+listak[7])
K[5,4]=-(4*listak[5]+4*listak[9]+listak[7])
K[5,6]=-(4*listak[5]+4*listak[9]+listak[7])
K[6,5]=-(4*listak[5]+4*listak[9]+listak[7])
K[6,7]=-(4*listak[5]+4*listak[9]+listak[7])
K[7,6]=-(4*listak[5]+4*listak[9]+listak[7])
K[7,7]=4*listak[5]+4*listak[9]+2*listak[7]+4*listak[8]
K[7,8]=-(4*listak[8]+listak[7])
K[8,7]=-(4*listak[8]+listak[7])
i+=1
while i<=10:
    K[i,i]=8*listak[8]+2*listak[7]
    K[i,i+1]=-(4*listak[8]+listak[7])
    K[i+1,i]=-(4*listak[8]+listak[7])
    i+=1
K[i,i]=4*listak[8]+5*listak[7]+listak[6]
K[i,i+1]=-(4*listak[7]+listak[6])
K[i+1,i]=-(4*listak[7]+listak[6])
i+=1
while i<=14:
    K[i,i]=8*listak[7]+2*listak[6]
    K[i,i+1]=-(4*listak[7]+listak[6])
    K[i+1,i]=-(4*listak[7]+listak[6])
    i+=1
K[i,i]=4*listak[7]+5*listak[5]+listak[6]
K[i,i+1]=-(5*listak[5])
K[i+1,i]=-(5*listak[5])
i+=1
while i<=19:
    K[i,i]=10*listak[5]
    i+=1
K[0,1]=-(6*listak[5]+2*listak[6]+4*listak[9]+listak[8])
K[1,0]=-(6*listak[5]+2*listak[6]+4*listak[9]+listak[8])
i=14
while i<=18:
    K[i,i+1]=-(5*listak[5])
    K[i+1,i]=-(5*listak[5])
    i+=1

#Matriz de masa expresada en toneladas (ton)
i = 0
mass = sp.zeros((20,20))
for piso in mass:
    if i < 4:
        mass[i,i]=6.2*3.5*12
    
    elif i>=4 and i<=7:
        mass[i,i]=6.2*3.5*8
    elif i>=8 and i <=19:
        mass[i,i]=6.2*3.5*4
    i+=1



f1 = 0.2 #frecuencia en hz
f2 = 2 #frecuencia en hz
chi1 = 0.025 #Coeficiente Rayleigh 2,5%
chi2 = 0.025 #COoficiente de Rayleigh 2,5%
a0 = 4*math.pi*f1*f2*(f1*chi2-f2*chi1)/(f1**2-f2**2)
a1 = (f1*chi1-f2*chi2)/(math.pi*(f1**2-f2**2))

#Matriz de Coeficientes expresado en Ggr/s , 'Giga gramos por segundo'
C = a0*mass + a1*K
identidad=np.identity(20) #Esquina superior izquierda de la matriz, Matriz identidad, 20x20
zeros=np.zeros((20,20)) #ESquina superior derecha de la matriz , MAtriz de ceros 20x20
invmass=inv(mass)#Se invierte la matriz de masas para los calculos siguientes
esqinfi=-(invmass.dot(K)) #Esquina inferior izquierda de la matriz M-1*K, 20x20
esqinfd=-(invmass.dot(C)) #ESquina inferior izquierda de la matriz M-1*C, 20x20

#Luego vamos a conglomerar las 4 matricez para generar la matriz A, 40x40
a1=np.concatenate((zeros,identidad), axis=1)
a2=np.concatenate((esqinfi,esqinfd), axis=1)

A=np.concatenate((a1,a2), axis=0)

sp.savez('mck.npz',M=mass,C=C,K=K,A=A) #Se guardan las matricez en una carpeta, en formato .npz
    