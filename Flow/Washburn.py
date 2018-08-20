# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 17:48:21 2018

@author: andrew
"""

import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
import pandas as pd

class Fluid:
    def __init__(self,viscosity,surface_tension):
        self.viscosity = viscosity
        self.surface_tension = surface_tension

# http://www.ddbst.com/en/EED/PCP/VIS_C516.php
hex_viscosity = 3 ## mPaS
hex_surface_tension = 27.5 ## mN/m


hexadecane = Fluid(hex_viscosity*1e-3, hex_surface_tension*1e-3)

def ellipse_circ1(a,b):  ##Ramanujan
    h = (a-b)**2 / (a+b)**2
    return pi*(a+b)*(1 + 3*h/(10 + np.sqrt(4 - 3*h)))

def washburn_eqn_Diameter(slope,gamma,eta,contact_angle):
    D = 4*eta*slope/(gamma*np.cos(contact_angle))
    return D

def Hagen_volume_flow_ellipse(dp,L,a,b,viscosity): # a,b semi axis
    x = np.pi * dp * a**3 * b**3
    y = 4*viscosity*(a**2 + b**2)
    return x/y

def Washburn_ellipse_slope(a,b,viscosity,contact_angle,surface_tension):
    c = ellipse_circ1(a,b)
    x = a * b * surface_tension * c * np.cos(contact_angle)
    y = 2 * viscosity * (a**2 + b**2)
    return x/y

def Washburn_flow_ellipse(ts,a,b,fluid,contact_angle):
    viscosity = fluid.viscosity
    surface_tension = fluid.surface_tension
    
    c = ellipse_circ1(a,b)
    
    L_sq = 1/(2*pi)*a*b*c/(a**2 + b**2)*surface_tension*np.cos(contact_angle)/viscosity*ts
    
    return np.sqrt(L_sq)

def Washburn_flow_circle(ts,r,fluid,contact_angle):
    viscosity = fluid.viscosity
    surface_tension = fluid.surface_tension
    
    L_sq = surface_tension*r*np.cos(contact_angle)*ts/(2*viscosity)
    
    return np.sqrt(L_sq)

def Washburn_const_ellipse(a,b,fluid,contact_angle):
    viscosity = fluid.viscosity
    surface_tension = fluid.surface_tension
    
    c = ellipse_circ1(a,b)
    
    W = 1/(2*pi)*a*b*c/(a**2 + b**2)*surface_tension*np.cos(contact_angle)/viscosity
    return W

def Washburn_flow_dist(a,b,viscosity,contact_angle,surface_tension,t):
    c = ellipse_circ1(a,b)
    x = a * b * surface_tension * c * np.cos(contact_angle) * t
    y = 2 * viscosity * (a**2 + b**2)
    
def surface_tension(perimeter,surface_tension, contact_angle):
    """ surface_tension(perimeter, surface_tension, contact_angle)
    Surface teension provided by constant cross-section. """
    
    return perimeter*surface_tension*contact_angle