# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 10:25:35 2016

@author: juanma
"""

#%%
# Librerías numéricas
from numpy import *
from matplotlib.pyplot import *
from matplotlib.mlab import find

#%% Función segmenta
def segmenta(Vector, Condicion = 0.5, Modo = "m"):
    '''
    segmenta
    ========
    
       Hace una segmentación inteligente para permitir
       la clasificacion de los puntos
       
    uso:
    ----
       
       >>> y = segmenta(Vector, Condicion, Modo)
    
           * *Vector*  : Vector de muestras equiespaciadas en el tiempo cuyos puntos queremos caracterizar
           
           * *Valor*   : El valor que determina la calidad del ajuste para el calculo de los maximos y minimos
           
           * *Modo*    : El procedimiento para calcular la calidad del ajuste. Puede tener dos valores
           
                     - 'i': Con lo que Valor indica el numero de iteraciones que queremos que haga el sistema
                     
                     - 'e': Con lo que Valor indica el error cuadratico medio maximo asumible.
                     
                     - 'm': Valor indica el maximo error asumible respecto la curva.

           *  y : Matriz con dos columnas y tantas filas como puntos de interes.
                  las columnas estan ordenadas por pares (x,y).
    
    :Version: 1.2.
       
    :Date: Noviembre de 2002
       
    :Author: Juan Manuel Belda Lois para el IBV.

    :Contact: jmbeldalois "at" gmail.com
    
    :Organization: IBV (Instituto de Biomecánica de Valencia).
    
    Traducido a python el 27 de Septiembre de 2016
    '''

    sgm = [] # Los segmentos
    
    # Hacemos un cast de Vector como array
    Vector = array(Vector)
    
    # La longitud del Vector
    lvect = len(Vector)

    # Ponemos los extremos del vector como primeros puntos
    sgm.append([0, Vector[0]])
    sgm.append([lvect-1, Vector[-1]])
    
    # Amplitud de la señal
    maxi = max(Vector)
    mini = min(Vector)
    amplitud = maxi - mini

    # Calculamos los errores
    asgm = array(sgm)
    asgm = interp(range(lvect),asgm[:,0],asgm[:,1])
    el_error = (Vector - asgm)
    Distancias = abs(el_error)
    
    # Calculamos el limite
    if Modo == "e":
        Limite = sum(el_error**2)/lvect
    elif Modo == "i":
        Limite = 2 * Condicion
    else:
        Limite = max(Distancias)/amplitud


    # Iteramos hasta que se cumple la condicion
    #Limite = Condicion + 1
    while(Limite > Condicion):
        # Buscamos la distancia máxima y añadimos el punto a la recta
        id_max = Distancias == max(Distancias)
        sgm.append([float(array(range(lvect))[id_max]), float(Vector[id_max])])
        
        # Ordenamos la salida
        sgm = sorted(sgm, key=lambda a_entry: a_entry[0]) 

        # Calculamos la distancia máxima y añadimos el punto a la recta
        asgm = array(sgm)
        vals = interp(range(lvect),asgm[:,0],asgm[:,1])
        el_error = Vector - vals
        Distancias = abs(el_error)
        
       
        # Calculamos el limite
        if Modo == "e":
            Limite = sum(el_error**2)/lvect
        elif Modo == "i":
            Limite -= 1
        else:
            Limite = max(Distancias)/amplitud
        
    
    return sgm


def clasifica_puntos(Vector, Valor=0.5, Modo="m", draw=False):
    '''
    clasifica_puntos
    ================
    
    Devuelve máximos y mínimos relativos, asi como los pasos por cero
    
    usos:
    -----
    
       >>> M, m, z = clasifica_puntos(Vector, Valor, Modo)
    
           * *Vector*  : Vector de muestras equiespaciadas en el tiempo cuyos puntos queremos caracterizar
           
           * *Valor*   : El valor que determina la calidad del ajuste para el calculo de los maximos y minimos
           
           * *Modo*    : El procedimiento para calcular la calidad del ajuste. Puede tener dos valores
           
                     - 'i': Con lo que Valor indica el numero de iteraciones que queremos que haga el sistema
                     
                     - 'e': Con lo que Valor indica el error cuadratico medio maximo asumible.
                     
                     - 'm': Valor indica el maximo error asumible respecto la curva.
    
           * *M* : Vector con los puntos que son maximos relativos
           
           * *m* : Vector con los puntos que son minimos relativos
           
           * *z*: Vector con los puntos que son paso por cero
    
           Admite un parámetro lógico (draw) para dibujar el vector y los
           puntos encontrados
    
    Ejemplos:
    ---------
    >>> # Calcula los máximos y mínimos relativos de vector
    >>> # tras 19 iteraciones.
    >>> maxi, mini, zeros = clasifica_puntos(vector,19,'i');
    >>> #                
    >>> # Calcula tras ajustar la curva para un error 
    >>> # cuadrático medio menor de 0.0018
    >>> # Adicionalmente dibuja la curva con los puntos característicos
    >>> maxi, mini, zeros = clasifica_puntos(total,18,'e',draw=True)
                      
    :Version: 1.3
    
    :Date: Diciembre de 2002.
    
    :Author: Juan Manuel Belda Lois
    
    :Contact: jmbeldalois "at" gmail.com
    
    :Organization: IBV (Instituto de Biomecánica de Valencia).
    
    Traducido a python el 27 de Septiembre de 2016
    '''
    
    puntos_interes = segmenta(Vector, Valor, Modo)
    
    l_vector = len(Vector) # La longitud del vector
    l_puntos = len(puntos_interes)
    
    ptos = array(puntos_interes) # Nos facilita la vida
    
    # Busqueda de máximos
    maximos = []
    minimos = []
    pasocero = []
    for c in arange(1, l_puntos-1, 1):
        ini = c - 1
        fin = c + 1
        
        pint = Vector[ptos[ini,0]:ptos[fin,0]]
        #pint = ptos[ini:fin,:] # El segment de interés actual
        
        # Vamos a ver si es un maximo
        if (ptos[c,1] >= ptos[ini,1]) & (ptos[c,1]>=ptos[fin,1]):
            imax = find(pint == max(pint))[0]
            #imax = find(pint[:,1] == max(pint[:,1]))
            maximos.append([ptos[ini,0]+imax, pint[imax]])
            #maximos.append(list(pint[imax[0],:]))
        # Vemos si es un mínimo
        elif (ptos[c,1] <= ptos[ini,1]) & (ptos[c,1]<=ptos[fin,1]):
            imin = find(pint == min(pint))[0]
            #imin = find(pint[:,1] == min(pint[:,1]))
            minimos.append([ptos[ini,0]+imin, pint[imin]])
            #minimos.append(list(pint[imin[0],:]))
        elif (sign(ptos[ini,1])!=sign(ptos[fin,1])):
            vals = abs(pint)
            #vals = abs(pint[:,1])
            iz = find(vals == min(vals))[0]
            pasocero.append([ptos[ini,0]+iz, pint[iz]])
            #pasocero.append(list(pint[iz[0],:]))

    amaxi = array(maximos)
    amini = array(minimos)
    apz = array(pasocero)
            
    if (draw):
        
        figure()
        plot(Vector)        
        if len(amaxi) > 0:
            plot(amaxi[:,0], amaxi[:,1],"^")
        if len(amini) > 0:
            plot(amini[:,0], amini[:,1],"v")
        if len(apz) > 0:
            plot(apz[:,0], apz[:,1],"d")
            
        show()
            
    return (amaxi, amini, apz)
        