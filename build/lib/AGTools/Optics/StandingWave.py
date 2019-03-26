# -*- coding: utf-8 -*-
"""
Created on Sun Oct  7 19:42:00 2018

@author: andrew
"""

import AGTools.Optics.BraggStackModel as BSM
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
pi = np.pi

ps = BSM.SellmeierMedium(BSM.ps_params)
si = BSM.SellmeierMedium(BSM.si_params)


def standing_wave_eqn(D, n1, n2, n3, w, z) :
    """
    D is thickness
    r12, r23 are reflectivity values between boundaries
    t12 is transmittivity at boundary
    tD is transmission of film
    k is propagation constant
    w is wavelength
    n is refractive index
    z is position in film
    """
    r12 = (n1 - n2)/(n1 + n2)
    r23 = (n2 - n3)/(n2 + n3)
    t12 = 2*n1/(n1 + n2)
    k2  = 2*pi*n2/w
    tD  = np.exp(-1j*k2*D)
    
    osc = np.exp(-1j*k2*z) + r23*tD**2*np.exp(1j*k2*z)
    ref = t12/(1 + r12*r23*tD**2)
    
    return osc*ref

def standing_wave_eqn_direct(D, t12, r12, r23, n2, w, z) :
    """
    D is thickness
    r12, r23 are reflectivity values between boundaries
    t12 is transmittivity at boundary
    tD is transmission of film
    k is propagation constant
    w is wavelength
    n is refractive index
    z is position in film
    """
    
    k2 = 2*pi*n2/w
    tD  = np.exp(-1j*k2*D)
    
    osc = np.exp(-1j*k2*z) + r23*tD**2*np.exp(1j*k2*z)
    ref = t12/(1 + r12*r23*tD**2)
    
    return osc*ref


def thickness_standing_wave_plot(w=405, thickness=1000, height=1000, ax=None):
    if not ax:
        ax = plt.subplot(111)
        
    film_thicknesses = np.arange(0,thickness)
    film_crosssection = [np.arange(0,i) for i in film_thicknesses]
    sws = [0.5*np.absolute(standing_wave_eqn(h, 1, ps.n(w), 1j, w, film_crosssection[h]))**2 for h in range(thickness)]
    sws_pad = [np.pad(a, (height-len(a),0), 'constant', constant_values=np.nan) for a in sws]
    sws_arr = np.array(sws_pad, dtype=float)
    
    max_amps = np.array([np.nanmax(a) for a in sws_arr])
    min_amps = np.array([np.nanmin(a) for a in sws_arr])
    
    ax.plot(film_thicknesses, max_amps)
    ax.plot(film_thicknesses, min_amps)
    ax.set_xlabel("Film Thickness (nm)")
    ax.set_ylabel("Peak Intensity")
    return ax

def thickness_standing_wave_heatmap_plot(w=405, thickness=1000, height=1000):
    

    film_thicknesses = np.arange(0,thickness)
    film_crosssection = [np.arange(0,i) for i in film_thicknesses]
    sws = [np.absolute(standing_wave_eqn(h, 1, ps.n(w), si.n(w), w, film_crosssection[h]))**2 for h in range(thickness)]
    sws_pad = [np.pad(a, (height-len(a),0), 'constant', constant_values=np.nan) for a in sws]
    sws_arr = np.array(sws_pad, dtype=float)
    
    max_amps = np.array([np.nanmax(a) for a in sws_arr])
    min_amps = np.array([np.nanmin(a) for a in sws_arr])
    
    cmap = matplotlib.cm.viridis
    cmap.set_bad('white',1.0)
    ax1 = plt.imshow(sws_arr.T, origin='upper', cmap=cmap, extent=(0,thickness,0,height))
    plt.colorbar()
    
    plt.xlabel("Film Thickness (nm)")
    plt.ylabel("Height (nm)")
    
    ax2 = plt.axes([0.28, 0.65, 0.2, 0.2])
    plt.xlabel("Film Thickness (nm)")
    plt.ylabel("Peak Intensity")
    
    ax2.plot(film_thicknesses, max_amps)
    #ax2.plot(film_thicknesses, min_amps/2)
    
    return sws_arr