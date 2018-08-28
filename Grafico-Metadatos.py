#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 19:34:00 2018

@author: jespildora
"""
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import glob
path = "/home/jespildora/Documents/Modelos Computacionales en IOC/sismos 0,2-0,8/*.txt"
files = glob.glob(path) 
for name in files:
        a=open(name)
        l1 = a.readline()
        l2= a.readline()
        l3= a.readline()
        l4= a.readline()
        l5= a.readline()
        l6= a.readline() 
        a = sp.loadtxt(name)*50
         
        dt = 1./200.
        Nt = a.size
         
         
        t = sp.arange(0, dt*Nt, dt)
         
        Ia = sp.zeros(Nt)
        v = sp.zeros(Nt)
        d = sp.zeros(Nt)
         
        v[1:] = sp.cumsum(a[1:] + a[0:-1])*dt/2
        d[1:] = sp.cumsum(v[1:] + v[0:-1])*dt/2
         
         
         
        g = 9.806
         
        a2 = a**2
        da2 = (a2[0:-1] + a2[1:])*dt/2
         
        Ia[1:] = sp.cumsum(da2)*sp.pi/(2*g)
         
        Ia_inf = Ia.max()
         
        i_PGA = sp.argmax(abs(a))
        t_PGA = t[i_PGA]
        PGA = (a[i_PGA])
         
        i_PGV = sp.argmax(abs(v))
        t_PGV = t[i_PGV]
        PGV = (v[i_PGV])
         
        i_PGD = sp.argmax(abs(d))
        t_PGD = t[i_PGD]
        PGD = (d[i_PGD])
         
         
        i_05 = sp.argmin( abs(Ia - 0.05*Ia_inf) )
        i_95 = sp.argmin( abs(Ia - 0.95*Ia_inf) )
         
        t_05 = t[i_05]
        Ia_05 = Ia[i_05]
         
        t_95 = t[i_95]
        Ia_95 = Ia[i_95]
         
        D_5_95 = t_95 - t_05
         
        
        metadatos = {}
        metadatos['fecha']=l1[20:30]
        metadatos['hora']=l1[31:40]
        metadatos['Lat']=l5[10:18]
        metadatos['Long']=l5[30:36]
        metadatos['Estacion nombre']=l4[11:16]
        metadatos['PGA']=PGA/10
        metadatos['PGV']=PGV*100
        metadatos['PGD']=PGD*100
        metadatos['Duracion']=D_5_95
        
         
        plt.figure().set_size_inches([9,6])
         
        plt.subplot(3,1,1)
        plt.plot(t,a/g)
        plt.axvline(t_05, color="k", linestyle="--")
        plt.axvline(t_95, color="k", linestyle="--")
        plt.text(t_PGA,PGA/g, "PGA={0:0.3f}g".format(abs(PGA)/g))
        plt.plot(t_PGA,PGA/g, "ob")
        plt.ylim([-PGA/8.,PGA/8.])
        plt.grid(True)
        plt.ylabel("Acc, $a$ (g)")
         
        plt.subplot(3,1,2)
        plt.plot(t,v*100)
        plt.axvline(t_05, color="k", linestyle="--")
        plt.axvline(t_95, color="k", linestyle="--")
        plt.text(t_PGV,PGV*100,"PGV={0:0.3f}cm/s".format(abs(PGV)*100))
        plt.plot(t_PGV,PGV*100, 'ob')
        plt.ylim([-PGV*150,PGV*150])
        plt.grid(True)
        plt.ylabel("Vel, $v$ (cm/s)")
         
         
        plt.subplot(3,1,3)
        plt.plot(t,d*100)
        plt.axvline(t_05, color="k", linestyle="--")
        plt.axvline(t_95, color="k", linestyle="--")
        plt.text(t_PGD,PGD*100, "PGD={0:0.3f}cm".format(abs(PGD)*100))
        plt.plot(t_PGD,PGD*100, "ob")
        plt.ylim([-PGD*150,PGD*150])
        plt.grid(True)
        plt.xlabel("Tiempo, $t$ (s)")
        plt.ylabel("Dis, $d$ (cm)")
         
        plt.subplot(3,1,1)
        plt.title('Fecha:  ' + metadatos['fecha'] + '  Hora:  ' + metadatos['hora'] + '  Estacion:  ' + metadatos['Estacion nombre'] + "    $D_{{5-95}} = {0:5.2f}$s".format(D_5_95))
         
        

plt.show()
