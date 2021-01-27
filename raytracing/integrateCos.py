#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 27 08:59:01 2021

@author: raphael
"""

import numpy as np

def integrateCos():
    N = 10000
    sigma = 0.25
    s = 0
    for i in range(N):
        u1 = np.random.rand()
        u2 = np.random.rand()
        xi = sigma*np.cos(2*np.pi*u1)*np.sqrt(-2*np.log(u2))
        p = 1/(sigma*np.sqrt(2*np.pi)) * np.exp(-xi*xi / (2*sigma*sigma)) # Loi de probabilité gaussienne
        s += np.cos(xi)**10/p/N # x_i
    return s

print(integrateCos())