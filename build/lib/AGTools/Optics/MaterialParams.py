# -*- coding: utf-8 -*-
"""
Created on Mon Aug 20 18:38:05 2018

@author: andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os.path as path


## Params: Sellmeier A, Sellmeier B, Sellmeier C, d_h, dummy, porosity, n_si,
##         dummy(window), d_l, n_air
## filmetrics.com / Silicon

ps_params = {'A': 1, 'B': 1.4435, 'C': 20216}
si_params = {'A': 1, 'B': 10.668, 'C': 90912}
n_air = 1
n_si = 3.67
df_si = pd.read_csv(path.join(path.split(__file__)[0],"Si.txt"),header=0,sep='\t')