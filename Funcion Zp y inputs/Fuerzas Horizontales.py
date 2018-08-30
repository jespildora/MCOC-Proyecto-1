#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 18:59:22 2018

@author: jespildora
"""
import scipy as sp
from numpy.linalg import inv

datos=sp.load('mck.npz')
C=datos['C']
M=datos['M']
K=datos['K']
A=datos['A']
invM=inv(M)
masa=sp.zeros(20)
datos = sp.load("registro_02.npz")
ac=datos['a']
p=0
while p<=19:
    masa[p]=M[p,p]
    p+=1

Mtotal=sp.sum(masa)
vs=max(ac)*Mtotal
Htotal=57.2
T=0.047*57.2**0.9
k=0.75+0.5*T
cv=sp.zeros(20)
j=1
hi=sp.zeros(20)
hi[0]=4.0
delta=sp.zeros(20)
delta[0]=0.04
while j<=19:
    hi[j]=(2.8+hi[j-1])
    delta[j]=0.001*hi[j]
    j+=1
j=0

hik=hi**k
mihi=sp.dot(masa,hik)


j=0

while j<=19:    
    cv[j]=(masa[j]*hik[j])/(mihi)    
    j+=1
Fi=sp.zeros(20)
Fi=cv*vs
Fdesp=K*delta
discipadores=sp.zeros(20)
discipadores=abs(Fi-Fdesp)/10
