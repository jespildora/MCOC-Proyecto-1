import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt

datos = sp.loadtxt("20161009-010221-C18O-HNE.txt")
x = sp.arange(len(datos))
f = interpolate.interp1d(x, datos)
xnew = sp.arange(0, len(datos), 1)
ynew = f(xnew)
plt.plot(x,datos,"o",xnew,ynew)
plt.show()