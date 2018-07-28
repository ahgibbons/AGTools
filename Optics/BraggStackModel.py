# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 14:14:52 2018

@author: andrew
"""

"""Multi Layer Bragg Model"""

import matplotlib.pyplot as plt
import numpy as np
exp = np.exp
sqrt = np.sqrt
from math import pi as PI
from numpy.linalg import det,matrix_power
from scipy.optimize import curve_fit
import pandas as pd

## Params: Sellmeier A, Sellmeier B, Sellmeier C, d_h, dummy, porosity, n_si,
##         dummy(window), d_l, n_air

ps_params = {'A': 1, 'B': 1.4435, 'C': 20216}
si_params = {'A': 1, 'B': 10.668, 'C': 90912}
n_air = 1
n_si = 3.67
df_si = pd.read_csv("Si.txt",header=0,sep='\t')

## filmetrics.com / Silicon


params = [1,1.4435,20216,100,0,0.25,3.67,1,100,1]
wavelengths = np.linspace(200,1000,1000)

class SellmeierLayer:
    def __init__(self, mat_params, thickness):
        self.A = mat_params['A']
        self.B = mat_params['B']
        self.C = mat_params['C']
        self.thickness = thickness
        
        
    def matrix(self,w):
        n = np.sqrt(self.A + ((self.B * w**2) / (w**2 - self.C)))
        
        k = 2*PI*n / w
        theta = k*self.thickness
        
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
        return M
    
class DataFrameLayer:
    def __init__(self, material_dataframe, thickness):
        self.dataframe = material_dataframe
        self.thickness = thickness
        
    def matrix(self,w):
        n = np.interp(w,self.dataframe['Wavelength(nm)'],self.dataframe['n'])
        k = 2*PI*n / w
        theta = k*self.thickness
        
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
    
class PorousLayer:
    def __init__(self, mat_params, porosity, thickness):
        self.A = mat_params['A']
        self.B = mat_params['B']
        self.C = mat_params['C']
        self.porosity = porosity
        self.thickness = thickness

    def matrix(self,w):
        nh = np.sqrt(self.A + ((self.B * w**2) / (w**2 - self.C)))
        n = 1 + ((nh -1 ) * (1-self.porosity)) 
        
        k = 2*PI*n / w
        theta = k*self.thickness
        
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
        return M

class SimpleLayer:
    def __init__(self, refractive_index, thickness):
        self.refractive_index = refractive_index
        self.thickness = thickness
        
    def matrix(self,w):
        k = 2*PI*self.refractive_index / w
        theta = k*self.thickness
        
        n = self.refractive_index
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
        return M
    
def curve_diff(ys1,ys2): ## Two curves of equal length with same xs domain
    diff = ys1 - ys2
    diffsq = diff * diff
    return sum(diffsq)
    
def gen_multilayer_matrix(w,layers,n_Left,right_params):
    BLinv = 0.5*np.array([[1,-1/n_Left],[1,1/n_Left]])
    A = right_params['A']
    B = right_params['B']
    C = right_params['C']
    n_si_sell = np.sqrt(A + ((B * w**2) / (w**2 - C)))
    BR = np.array([[1,1],[-1*n_si_sell,n_si_sell]])
    
    MT = BLinv
    
    for layer in layers:
        m = layer.matrix(w)
        MT = MT @ m
    
    MT = MT @ BR
    
    r = MT[1][0] / MT[0][0]
    
    return (r * r.conj()).real

def gen_multilayer_matrix_df_substrate(w,layers,n_Left,df_right):
    n_Right = np.interp(w,df_right['Wavelength(nm)'],df_right['n'])
    
    BLinv = 0.5*np.array([[1,-1/n_Left],[1,1/n_Left]])
    BR = np.array([[1,1],[-1*n_Right,n_Right]])
    
    MT = BLinv
    
    for layer in layers:
        m = layer.matrix(w)
        MT = MT @ m
    
    MT = MT @ BR
    
    r = MT[1][0] / MT[0][0]
    
    return (r * r.conj()).real

def gen_multilayer_matrix_fixed_substrate(w,layers,n_Left,n_Right):
    BLinv = 0.5*np.array([[1,-1/n_Left],[1,1/n_Left]])
    BR = np.array([[1,1],[-1*n_Right,n_Right]])
    
    MT = BLinv
    
    for layer in layers:
        m = layer.matrix(w)
        MT = MT @ m
    
    MT = MT @ BR
    
    r = MT[1][0] / MT[0][0]
    
    return (r * r.conj()).real
    
def cos_layer(ws,mat_params,porosity,d1,d2,nlayers,n_Left,n_Right):
    layerSolid = SellmeierLayer(mat_params,d1)
    layerPorous = PorousLayer(mat_params,porosity,d2)
    layers = [layerSolid,layerPorous]*nlayers + [layerSolid]
    
    rs = np.array([gen_multilayer_matrix(w,layers,n_Left,n_Right) for w in ws])
    return rs

def cos_layer_fitter(ws,porosity,d1,d2):
    return cos_layer(ws,ps_params,porosity,d1,d2,5,n_air,n_si)

    
def multilayer_matrix(n1,n2,L1,L2,nlayers,w):
    B0inv = 0.5*np.array([[1,-1],[1,1]])
    Bsub  = np.array([[1,1],[-1*n_si,n_si]])
    k1 = 2*PI*n1 / w
    k2 = 2*PI*n2 / w
    
    theta1 = k1*L1
    theta2 = k2*L2
    
    
    M1 = np.array([[np.cos(theta1), 1j*np.sin(theta1)/n1],
                    [1j*n1*np.sin(theta1), np.cos(theta1)]])
    
    M2 = np.array([[np.cos(theta2), 1j*np.sin(theta2)/n2],
                    [1j*n2*np.sin(theta2), np.cos(theta2)]])
 
    
    M1M2 = np.matmul(M1,M2)
    
    MN = matrix_power(M1M2,nlayers)
    
    MT = B0inv @ MN @ M1 @ Bsub
    
    r = MT[1][0] / MT[0][0]

    return (r*r.conj()).real    


def porous_multilayer(material_params,layers,porosity,d_thick,
                      d_porous,n_medium,n_reflector,x): ## w is params vector, x is wavelength, n is number of layers
    
    A = material_params['A']
    B = material_params['B']
    C = material_params['C']
    
    jj = 0 + 1j
    n_h = np.sqrt(A + ((B * x**2) / (x**2 - C)))
    n_l = 1 + ((n_h -1 ) * porosity)

    phi_h  = (2*PI * n_h * d_thick)/x
    phi_l  = (2*PI * n_l * d_porous)/x
    r_lh   = (-n_h + n_l) / (n_h + n_l)
    r_hSi  = (n_h - n_reflector) / (n_h + n_reflector)
    r_hl   = (n_h - n_l) / (n_h+n_l)
    r_hAir = (n_h - n_medium) / (n_h + n_medium)
    r_Airh = (n_medium - n_h) / (n_h + n_medium)
    
    r_H    = (r_lh + r_hSi * exp(-jj*2*phi_h)) / (1 - r_hl * r_hSi*exp(-jj*2*phi_h)) ## Bottom Layer
    r_L    =  (r_hl + r_H * exp(-jj*2*phi_l)) / (1 - r_lh * r_H*exp(-jj*2*phi_l))    ## Bottom Layer
    
    for i in range(layers-1):
        r_H = (r_lh + r_L * exp(-jj*2*phi_h))/(1 - r_hl * r_L*exp(-jj*2*phi_h))
        r_L = (r_hl + r_H * exp(-jj*2*phi_l))/(1 - r_lh * r_H*exp(-jj*2*phi_l))
    
    r_H = (r_Airh + r_L * exp(-jj*2*phi_h))/(1 - r_hAir * r_L*exp(-jj*2*phi_h))
	
    r_Hconj = (r_H).conj()
    rr_H    = (r_H*r_Hconj).real

    return rr_H

## Params: Sellmeier A, Sellmeier B, Sellmeier C, d_h, dummy, porosity, n_si,
##         dummy(window), d_l, n_air
def multi_film1(w,x): ## w is params vector, x is wavelength
    
    jj = 0 + 1j
    n_h=sqrt(w[0]+((w[1]*x**2)/(x**2-w[2])))
    n_l=1+((n_h-1)*w[5])

    fai_h = (2*PI*n_h*w[3])/x
    fai_l =  (2*PI*n_l*w[8])/x
    r_lh = (-n_h+n_l)/(n_h+n_l)
    r_hSi =  (n_h-w[6])/(n_h+w[6])
    r_hl = (n_h-n_l)/(n_h+n_l)
    r_hAir = (n_h-w[9])/(n_h+w[9])
    r_Airh = (w[9]-n_h)/(n_h+w[9])
    r_H = (r_lh+r_hSi*exp(-jj*2*fai_h))/(1-r_hl*r_hSi*exp(-jj*2*fai_h)) ## Bottom Layer
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))    ## Bottom Layer
    r_H = (r_lh+r_L *exp(-jj*2*fai_h))/(1-r_hl*r_L *exp(-jj*2*fai_h))
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))
    r_H = (r_lh+r_L *exp(-jj*2*fai_h))/(1-r_hl*r_L *exp(-jj*2*fai_h))
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))
    r_H = (r_lh+r_L *exp(-jj*2*fai_h))/(1-r_hl*r_L *exp(-jj*2*fai_h))
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))
    r_H = (r_lh+r_L *exp(-jj*2*fai_h))/(1-r_hl*r_L *exp(-jj*2*fai_h))
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))
    r_H = (r_lh+r_L *exp(-jj*2*fai_h))/(1-r_hl*r_L *exp(-jj*2*fai_h))
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))
    r_H = (r_lh+r_L *exp(-jj*2*fai_h))/(1-r_hl*r_L *exp(-jj*2*fai_h))
    r_L =  (r_hl+r_H*exp(-jj*2*fai_l))/(1-r_lh*r_H*exp(-jj*2*fai_l))
    
    r_H = (r_Airh+r_L *exp(-jj*2*fai_h))/(1-r_hAir*r_L *exp(-jj*2*fai_h))
	
    r_Hconj = (r_H).conj()+w[4]
    rr_H = w[7]*(r_H*r_Hconj).real

    return rr_H

def compare(mparams,dh,dl,porosity,n_medium,n_substrate,ws):
    param_list = [mparams['A'],mparams['B'],mparams['C'],dh,0,porosity,
                  n_substrate,1,dl,n_medium]
    masateru = multi_film1(param_list,ws)
    layerps = SellmeierLayer(mparams,dh)
    layerpore = PorousLayer(mparams, porosity,dl)
    layers = [layerps,layerpore]*7 + [layerps]
    matrix = [gen_multilayer_matrix(layers,n_medium,n_substrate,w) for w in wavelengths]
    
    plt.plot(ws,masateru)
    plt.plot(ws,matrix)