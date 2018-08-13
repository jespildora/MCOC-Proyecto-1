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
e=23500.0
listai=[(1.0/3.0)*(60.0**4.),(1.0/3.0)*(70.0**4.0),(1.0/3.0)*(80.0**4.0),(1.0/3.0)*(90.0**4.0),(1.0/3.0)*(100.0**4.0)]
listak=[(12.0*e*listai[0])/4.0,(12.0*e*listai[1])/4.0,(12.0*e*listai[2])/4.0,(12.0*e*listai[3])/4.0,(12.0*e*listai[4])/4.0,(12.0*e*listai[0])/2.8,(12.0*e*listai[1])/2.8,(12.0*e*listai[2])/2.8,(12.0*e*listai[3])/2.8,(12.0*e*listai[4])/2.8]

K = np.zeros((20,20))
K[0,0]=6*listak[0]+2*listak[1]+4*listak[4]+listak[3]+6*listak[5]+2*listak[6]+4*listak[9]+listak[8]
i=1
while i<=2:
    K[i,i]=12*listak[5]+4*listak[6]+8*listak[9]+2*listak[8]
    i+=1
K[3,3]=10*listak[5]+2*listak[6]+8*listak[9]+listak[7]+listak[8]
while i<=6:
    K[i,i]=8*listak[5]+8*listak[9]+2*listak[7]
    i+=1
K[7,7]=4*listak[5]+4*listak[9]+2*listak[7]+4*listak[8]
while i<=10:
    K[i,i]=8*listak[8]+2*listak[7]
    i+=1
K[11,11]=4*listak[8]+5*listak[7]+listak[6]
while i<=14:
    K[i,i]=8*listak[7]+2*listak[6]
    i+=1
K[15,15]=4*listak[7]+5*listak[5]+listak[6]
while i<=19:
    K[i,i]=10*listak[5]
    i+=1
print K

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
print mass