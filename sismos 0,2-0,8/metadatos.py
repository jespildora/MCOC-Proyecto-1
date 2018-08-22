import scipy as sp
from scipy import interpolate
import matplotlib.pyplot as plt

datos = open("20140403-024228-T08A-HNE.txt")
l1 = datos.readline()
l2= datos.readline()
l3= datos.readline()
l4= datos.readline()
l5= datos.readline()
l6= datos.readline()
l7= datos.readline()
l8 = datos.readline()
l9= datos.readline()
l10 = datos.readline()
l11 = datos.readline()
l12 = datos.readline()

dic = {}
dic['fecha']=l1[20:30]
dic['hora']=l1[31:40]
dic['Epi_Lat']=l5[10:18]
dic['Epi_Long']=l5[30:36]
dic['Estacion nombre']=l4[11:16]
dic['Estacion_Lat']=l7[15:]
dic['Estacion_Long']=l8[15:]
dic['Profundidad']=l10[13:]
dic['Magnitud']=l11[4:]
dic['Componente']=l9[12:]
print dic