# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 18:38:05 2018

@author: andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os.path as path
from io import StringIO
import AGTools.Optics.PET as PET
from scipy.optimize import curve_fit

## Params: Sellmeier A, Sellmeier B, Sellmeier C, d_h, dummy, porosity, n_si,
##         dummy(window), d_l, n_air
## filmetrics.com / Silicon

fit_eqn = lambda x, b, c: 1 + b*x**2/(x**2 - c)

ps_params = {'A': 1, 'B': 1.4435, 'C': 20216}
si_params = {'A': 1, 'B': 10.668, 'C': 90912}
pc_params = {'A': 1, 'B': 1.4182, 'C': 21304} # https://refractiveindex.info/?shelf=organic&book=polycarbonate&page=Sultanova
# PET https://www.filmetrics.com/refractive-index-database/PET/Estar-Melinex-Mylar
df_pet = pd.read_csv(StringIO(PET.PET_data),header=0,sep='\t')
pet_fit = curve_fit(fit_eqn, df_pet['Wavelength(nm)'], df_pet['n'])

pet_params = {'A': 1, 'B': pet_fit[0][0], 'C': pet_fit[0][1]}


n_air = 1
n_si = 3.67
df_si = pd.read_csv(path.join(path.split(__file__)[0],"Si.txt"),header=0,sep='\t')