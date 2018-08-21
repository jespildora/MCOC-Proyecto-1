import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import glob
path = "/home/erenchun/Desktop/MCOC-Proyecto-1/sismos 0,2-0,8/*.txt"
files = glob.glob(path) 
for name in files:
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
    print name
    print "t_95 = ", t_95
    print "Ia_95 = ", Ia_95
    print "D_5_95 = ", D_5_95
 
    plt.figure().set_size_inches([9,6])
 
    plt.subplot(3,1,1)
    plt.plot(t,a/g)
    plt.axvline(t_05, color="k", linestyle="--")
    plt.axvline(t_95, color="k", linestyle="--")
    plt.text(t_PGA,PGA/g, "PGA={0:0.3f}g".format(abs(PGA)/g))
    plt.plot(t_PGA,PGA/g, "ob")
    plt.ylim([-0.6,0.6])
    plt.grid(True)
    plt.ylabel("Acc, $a$ (g)")
 
    plt.subplot(3,1,2)
    plt.plot(t,v*100)
    plt.axvline(t_05, color="k", linestyle="--")
    plt.axvline(t_95, color="k", linestyle="--")
    plt.text(t_PGV,PGV*100, "PGV={0:0.3f}cm/s".format(abs(PGV)*100))
    plt.plot(t_PGV,PGV*100, "ob")
    plt.ylim([-15,15])
    plt.grid(True)
    plt.ylabel("Vel, $v$ (cm/s)")
 
 
    plt.subplot(3,1,3)
    plt.plot(t,d*100)
    plt.axvline(t_05, color="k", linestyle="--")
    plt.axvline(t_95, color="k", linestyle="--")
    plt.text(t_PGD,PGD*100, "PGD={0:0.3f}cm".format(abs(PGD)*100))
    plt.plot(t_PGD,PGD*100, "ob")
    plt.ylim([-15,15])
    plt.grid(True)
    plt.xlabel("Tiempo, $t$ (s)")
    plt.ylabel("Dis, $d$ (cm)")
 
    plt.subplot(3,1,1)
    plt.title(name[:-4] + "    $D_{{5-95}} = {0:5.2f}$s".format(D_5_95))
 
    plt.tight_layout()
 
    plt.show()
