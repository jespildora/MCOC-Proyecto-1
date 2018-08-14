# MCOC-Proyecto-1
MCOC-Proyecto-1

Introducción
==============
En este proyecto lo que se espera es lograr programar y conceptualizar un problema de dinamica de estructuras, con el fin de verificar la incertidumbre que genera el calculo de la estructura, mediante metodos numericos.
EL lenguaje que se usara, sera principalmente python, junto con excel.

Matrices M, C, K, A
==============
El primer archivo de python es el codigo de como creamos nuestras matrices (M C K A), primero las obtuvimos mediante ecxel donde lo que hizo fue cargar los distintos valores del edificio y sus dimensiones y asi crear la matriz K donde fuimos sumando y restando la distintas rigideces por cada piso correspondiente, al revisar nuestros datos y compararlos con valores obtenidos en sap2000 creamos las matrices usado codigo de python ya que asi vamosa poder tener estas matrices guardadas en un tipo de archivo .npz el cual nos servira para poder seguir con esta entrega mas adelante sin hacer el codigo tan engorroso y lento 

Matriz M: Matriz de masa por cada piso
Matriz K: Matriz de rigidez de la estructura
Matriz C: Suma ponderada de las matrices M y K
Matriz A: Concatenación de matrices operadas
