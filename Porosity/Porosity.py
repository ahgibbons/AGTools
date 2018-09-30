# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 12:34:09 2018

@author: andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def ellipse_circ1(a,b):  ##Ramanujan
    h = (a-b)**2 / (a+b)**2
    return pi*(a+b)*(1 + 3*h/(10 + np.sqrt(4 - 3*h)))
    
def washburn_D_ellipse(a,b): # L^2 = Dt
    circ = ellipse_circ1(a,b)
    return (a*b)*circ / (2*(a**2 + b**2))

def predict_flow(sem_df,viscosity,surface_tension):
    pass

def washburn_film(sem_df):
    return washburn_D_ellipse(sem_df['Major'],sem_df['Minor'])
