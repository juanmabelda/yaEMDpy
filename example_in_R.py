# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 09:56:47 2017

@author: juanma
"""

from numpy import *
from emd import *

#%% This is a sinthetic signal from the R library example
ndata = 3000
tt2 = linspace(0,9,ndata)
xt2 = sin(pi*tt2) + sin(2*pi*tt2) + sin(6*pi*tt2) + 0.5*tt2
plot(tt2, xt2)

#%% Extraemos las componentes
imf, res = newimf(tt2, xt2, draw=True)

#%%
figure()
imf, res = newimf(tt2, res, draw=True)

#%%
imfs, res = extractimf(tt2, xt2, nIter=3, draw=True)