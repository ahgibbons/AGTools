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
import os.path as path

import AGTools.Spectrometer.parseSpectral as PS
from AGTools.Optics.MaterialParams import ps_params, si_params, n_air, n_si, df_si



params = [1,1.4435,20216,100,0,0.25,3.67,1,100,1]
wavelengths = np.linspace(200,1000,1000)

#spec_file = open(r"C:\Users\andrew\Dropbox\Spectrometer\Andrew_spectraldata\20170418_PS\samples.csv",'r')
#spec_data = PS.genSpecData(spec_file.read())
#ps_spec = spec_data[1]

class SellmeierMedium:
    def __init__(self, params):
        self.params = params
        self.A = params['A']
        self.B = params['B']
        self.C = params['C']     
        
    def n(self,w):
        return np.sqrt(self.A + ((self.B * w**2) / (w**2 - self.C)))


def sellmeierEqn(A,B,C,w):
    n = np.sqrt(A + ((B * w**2) / (w**2 - C)))
    return n

def cauchyEqn(A,B,C,ws):
    ws2 = ws/1000
    return A + B/np.power(ws2,2) + C/np.power(ws2,4)

class SellmeierLayer:
    def __init__(self, mat_params, thickness, name=None):
        """
        SellmeierLayer(mat_params, thickness, name=None)
        """
        self.A = mat_params['A']
        self.B = mat_params['B']
        self.C = mat_params['C']
        self.thickness = thickness
        self.name = name
        
        
    def matrix(self,w):
        """
        matrix(w)
        Generate transmission matrix for layer for wavelength w
        """
        n = np.sqrt(self.A + ((self.B * w**2) / (w**2 - self.C)))
        
        k = 2*PI*n / w
        theta = k*self.thickness
        
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
        return M
    
class DataFrameLayer:
    def __init__(self, material_dataframe, thickness, name=None):
        self.dataframe = material_dataframe
        self.thickness = thickness
        self.name = name
        
    def matrix(self,w):
        n_real = np.interp(w,self.dataframe['Wavelength(nm)'],self.dataframe['n'])
        n_imag = np.interp(w,self.dataframe['Wavelength(nm)'],self.dataframe['k'])
        n = n_real + 1j*n_imag
        k = 2*PI*n / w
        theta = k*self.thickness
        
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
    
        return M
    
class PorousLayer:
    def __init__(self, mat_params, porosity, thickness, name = None):
        self.A = mat_params['A']
        self.B = mat_params['B']
        self.C = mat_params['C']
        self.porosity = porosity
        self.thickness = thickness
        self.name = name

    def matrix(self,w):
        nh = np.sqrt(self.A + ((self.B * w**2) / (w**2 - self.C)))
        n = 1 + ((nh -1 ) * (1-self.porosity)) 
        
        k = 2*PI*n / w
        theta = k*self.thickness
        
        M = np.array([[np.cos(theta), 1j*np.sin(theta)/n],
                        [1j*n*np.sin(theta), np.cos(theta)]])
        return M

class SimpleLayer:
    """
    SimpleLayer(refractive_index, thickness, name = None)
    Simple layer with constant refractive index.
    """
    def __init__(self, refractive_index, thickness, name = None):
        """
        SimpleLayer(refractive_index, thickness, name = None)
        Simple layer with constant refractive index.
        """
        self.refractive_index = refractive_index
        self.thickness = thickness
        self.name = name

    def __repr__(self):
        return "<SimpleLayer: n = {:.03f},  t = {:.03f}>".format(self.refractive_index, self.thickness)
        
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

def thin_film(ws,mat_params,thickness,substrate_df):
    layer = [SellmeierLayer(mat_params,thickness)]
    rs = [gen_multilayer_matrix_df_substrate(w,layer,n_air,substrate_df) for
          w in ws]
    return rs

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

def fit_porous_film_thickness(ws,ys,mat_params,porosity,df_sub,
                              nlayers,t0,t1):
    def fit_func(xs,t0,t1):
       layer_ps = SellmeierLayer(mat_params,t0)
       layer_pore = PorousLayer(mat_params,porosity,t1)
       layers = [layer_ps,layer_pore]*nlayers + [layer_ps]
       rs = [gen_multilayer_matrix_df_substrate(x,layers,n_air,df_sub) for
             x in xs]
       return rs
   
    a = curve_fit(fit_func,ws,ys,[t0,t1])
    ys_fit = fit_func(ws,a[0][0],a[0][1])
    
    return (a, ys_fit)
       

def gen_multilayer_matrix_df_substrate(w,layers,n_Left,df_right):
    n_Right = np.interp(w,df_right['Wavelength(nm)'],df_right['n'])
    k_Right = np.interp(w,df_right['Wavelength(nm)'],df_right['k'])
    ncomp_Right = n_Right + 1j*k_Right
    
    BLinv = 0.5*np.array([[1,-1/n_Left],[1,1/n_Left]])
    BR = np.array([[1,1],[-1*ncomp_Right,ncomp_Right]])
    
    MT = BLinv
    
    for layer in layers:
        m = layer.matrix(w)
        MT = MT @ m
    
    MT = MT @ BR
    
    r = MT[1][0] / MT[0][0]
    
    return (r * r.conj()).real

def gen_multilayer_matrix_fixed_substrate_reflectance(w,layers,n_Left,n_Right):
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
    t = 1 / MT[0][0]

    R = (r * r.conj()).real
    T = (t * t.conj()).real
    
    return (R,T)
    
def Rayleigh_Scattering(w,d,n,distance,scatter_angle=0):
    I1 = (1 + (np.cos(scatter_angle))**2)/(2*distance**2)
    I2 = np.power(2*PI/w,4)*(n**2 -1)**2/(n**2+2)**2*(d/2)**6
    return I1*I2

def cos_layer(ws,mat_params,porosity,d1,d2,nlayers,n_Left,df_Right):
    layerSolid = SellmeierLayer(mat_params,d1)
    layerPorous = PorousLayer(mat_params,porosity,d2)
    layers = [layerSolid,layerPorous]*nlayers + [layerSolid]
    
    rs = np.array([gen_multilayer_matrix_df_substrate(w,layers,n_Left,df_Right) for w in ws])
    return rs

def cos_layer_fitter(ws,porosity,d1,d2):
    return cos_layer(ws,ps_params,porosity,d1,d2,5,n_air,n_si)

def cos_5layer_fitter(ws,t0,t1,t2,t3,t4,t5,t6,p1,p3,p5):
    layer0 = SellmeierLayer(ps_params,t0)
    layer2 = SellmeierLayer(ps_params,t2)
    layer4 = SellmeierLayer(ps_params,t4)
    layer6 = SellmeierLayer(ps_params,t6)
    player1 = PorousLayer(ps_params,p1,t1)
    player3 = PorousLayer(ps_params,p3,t3)
    player5 = PorousLayer(ps_params,p5,t5)
    
    layers = [layer0,player1,layer2,player3,layer4,player5,layer6]
    
    rs = [gen_multilayer_matrix_df_substrate(w,layers,n_air,df_si) for w in ws]
    return rs
    
def cos_4layer_fitter(ws,t0,t1,t2,t3,t4,t5,t6,t7,t8,p1,p3,p5,p7):
    layer0 = SellmeierLayer(ps_params,t0)
    layer2 = SellmeierLayer(ps_params,t2)
    layer4 = SellmeierLayer(ps_params,t4)
    layer6 = SellmeierLayer(ps_params,t6)
    layer8 = SellmeierLayer(ps_params,t7)
    player1 = PorousLayer(ps_params,p1,t1)
    player3 = PorousLayer(ps_params,p3,t3)
    player5 = PorousLayer(ps_params,p5,t5)
    player7 = PorousLayer(ps_params,p7,t7)
    
    layers = [layer0,player1,layer2,player3,layer4,player5,layer6,player7,layer8]
    
    rs = [gen_multilayer_matrix_df_substrate(w,layers,n_air,df_si) for w in ws]
    return rs
    
def bilayer_fitter(ws,t0,t1,p):
    layer0 = SellmeierLayer(ps_params,t0)
    layer1 = PorousLayer(ps_params,p,t1)
    layers = [layer0,layer1]*4 + [layer0]
    
    rs = [gen_multilayer_matrix_df_substrate(w,layers,n_air,df_si) for w in ws]
    return rs

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
