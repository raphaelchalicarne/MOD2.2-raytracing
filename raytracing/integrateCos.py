#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 08:59:01 2021

@author: raphael
"""

import numpy as np

def isInInterval(x):
    return (-np.pi/2 < x < np.pi/2)

def integrateCos():
    N = 100000
    sigma = 1
    somme = 0
    for i in range(N):
        u1 = np.random.rand()
        u2 = np.random.rand()
        u3 = np.random.rand()
        u4 = np.random.rand()
        
        xi = sigma*np.cos(2*np.pi*u1)*np.sqrt(-2*np.log(u2))
        yi = sigma*np.sin(2*np.pi*u1)*np.sqrt(-2*np.log(u2))
        zi = sigma*np.cos(2*np.pi*u3)*np.sqrt(-2*np.log(u4))
        wi = sigma*np.sin(2*np.pi*u3)*np.sqrt(-2*np.log(u4))
        
        if isInInterval(xi) & isInInterval(yi) & isInInterval(zi) & isInInterval(wi):
            
            pxi = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-xi*xi / (2*sigma*sigma)) # Loi de probabilité gaussienne
            pyi = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-yi*yi / (2*sigma*sigma))
            pzi = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-zi*zi / (2*sigma*sigma))
            pwi = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-wi*wi / (2*sigma*sigma))
            
            somme += np.cos(xi + yi + zi + wi)**2/(pxi * pyi * pzi * pwi)/N
    return somme

print(integrateCos())