# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 10:00:18 2017

@author: juanma
"""

from numpy import *
from clasifica_puntos import *
from scipy.interpolate import interp1d

from scipy.signal import hilbert

def extractimf(x_val, y_val, nIter=1, draw=False):
    '''Extract a number of imfs'''
    imfs = []
    res = y_val
    
    for c in range(nIter):
        imf, res = newimf(x_val, res)
        imfs.append(imf)
        
    if draw:
        for c in range(nIter):
            subplot(nIter+1,1,c+1)
            plot(x_val, imfs[c])
            
        subplot(nIter+1,1,nIter+1)
        plot(x_val, res)
        
    return imfs, res


def newimf(x_val, y_val, draw=False):
    
    the_error = (max(y_val) - min(y_val))*0.05/max(shape(y_val))
    
    maxi, mini, _ = clasifica_puntos(y_val, the_error, "e")
        
    # If the first value doesn't exist we added the next maximum as value
    sup, idx_sup, maxi = _do_the_spline(x_val, y_val, maxi, kind="maximum")
    inf, idx_inf, mini = _do_the_spline(x_val, y_val, mini, kind="minimum")

    res = (sup + inf)/2.
    imf = y_val - res

    if draw:
        subplot(2,1,1)
        plot(x_val, y_val)
        plot(x_val[idx_sup], maxi[:,1],"o")
        plot(x_val, sup)
    
        plot(x_val[idx_inf], mini[:,1],"o")
        plot(x_val, inf)
        plot(x_val, res)
                
        subplot(2,1,2)
        plot(x_val, imf)
        plot(x_val, res)
        legend(["imf", "res"])
    
    return imf, res


    # Same for the minimum

def _do_the_spline(x_val, y_val, points, kind="maximum"):
    '''Perform the spline fitting to the found points'''
    
    if kind == "maximum":
        first_value = max([y_val[0], points[0,1]])
        last_value = max([y_val[-1], points[-1,1]])
    else:
        first_value = min([y_val[0], points[0,1]])
        last_value = min([y_val[-1], points[-1,1]])

    # If the first value doesn't exist we added the next maximum as value
    if points[0,0] != 0:
        points = vstack((array([0, first_value]), points))
        
    # The same for the last value
    if points[-1,0] != x_val[-1]:
        points = vstack((points, array([len(x_val)-1, last_value])))
        
    indexes = array(points[:,0],dtype=int)
    f_sup = interp1d(x_val[indexes], points[:,1],kind="cubic")#,fill_value=x_val)
    
    return f_sup(x_val), indexes, points
    
#analytic_signal = hilbert(signal)
#amplitude_envelope = np.abs(analytic_signal)
#instantaneous_phase = np.unwrap(np.angle(analytic_signal))
#instantaneous_frequency = np.diff(instantaneous_phase) / (2.0*np.pi) * fs
