import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt
import glob
path = "/home/jespildora/Desktop/MCOC-Proyecto-1/sismos 0,2-0,8/*.txt"
files = glob.glob(path) 

i=1
for name in files:
   
    a=open(name)
    l1 = a.readline()
    l2= a.readline()
    l3= a.readline()
    l4= a.readline()
    l5= a.readline()
    l6= a.readline() 
    l7= a.readline()
    l8 = a.readline()
    l9= a.readline()
    l10 = a.readline()
    l11 = a.readline()
    l12 = a.readline()
    a = sp.loadtxt(name)

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
    metadatos['Estacion_Lat']=l7[15:]
    metadatos['Estacion_Long']=l8[15:]
    metadatos['Profundidad']=l10[13:]
    metadatos['Magnitud']=l11[4:]
    metadatos['Componente']=l9[12:]
    metadatos['PGA']=PGA/10
    metadatos['PGV']=PGV*100
    metadatos['PGD']=PGD*100
    metadatos['Duracion']=D_5_95

    if i<=9:
        b="registro_0"+str(i)
        sp.savez(b,a=a,metadatos=metadatos,t=t)
        
    if i>=10:
        b="registro_"+str(i)
        sp.savez(b,a=a,metadatos=metadatos,t=t)
    print b
    print abs(metadatos['PGA'])
    i+=1
    CAP=[0, 0, 0, 0, 0, 500, 500, 0, 150, 150, 0, 0, 0, 0, 0, 500, 800, 800, 800, 800]