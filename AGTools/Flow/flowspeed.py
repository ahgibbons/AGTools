# -*- coding: utf-8 -*-
"""
Created on Sat May 19 20:07:50 2018

@author: andrew
"""

""" Channel flow calculations """

import matplotlib.pyplot as plt

import Sivaniah.Plot as Plot
from Sivaniah.SavitzkyGolay import savitzky_golay
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import numpy as np
import os
from PIL import Image
from PIL import ImageFilter

import numpy as np
from scipy.optimize import curve_fit

fit_func = (lambda t,a,n: a*t**n)

fit_sqrt = (lambda t,a: a*np.power(t,-0.5))


def calcAllPositions(imgseq,channelwidth,threshold=190,gaussian_radius=5,scale=1):
    positions = [calcPosition(im,channelwidth,threshold,gaussian_radius,scale) for 
                 im in imgseq]
    return positions

def calcPosition(imgpath, channelwidth, threshold=190,
                 gaussian_radius=5, scale=1, channel=2, fmt='HSV'):
    img = Image.open(imgpath)
    blurred = img.filter(ImageFilter.GaussianBlur(radius=gaussian_radius)).convert(fmt)
    v = blurred.split()[channel]
    threshold_img = v.point(lambda p: p > threshold and 255)
    area = threshold_img.histogram()[0]
    length = area/channelwidth
    return length*scale
    
def filteredGradient(ps,frame_interval=1,fps=0.5):
    time_interval = frame_interval/fps
    gs = np.gradient(ps,time_interval)
    ts = [i/fps*frame_interval for i in range(len(ps))]
    tg = zip(ts,gs)
    tg2 = [p for p in tg if p[1]>0]
    t,g = zip(*tg2)
    return (t,g)

@Plot.plotDecorator
def plot_with_fit(ax,ts,gs,markersize=50,markerstyle='o',
                  scattercolor='b',linecolor='r'):
    fits1 = curve_fit(fit_func,ts[1:],gs[1:])
    fits2 = curve_fit(fit_sqrt,ts[1:],gs[1:])
    
    xs_fit = [fit_func(t,fits1[0][0],fits1[0][1]) for t in ts[1:]]
    xs_sqrt = [fit_sqrt(t,fits2[0][0]) for t in ts[1:]]
    
    ax.scatter(ts,gs,color=scattercolor,s=markersize,marker=markerstyle,label="Hexadecane")
    #ax.plot(ts[1:],xs_fit,color='k',label="$at^{-n}$")
    ax.plot(ts[1:],xs_sqrt,color=linecolor,label="$bt^{-1/2}$")
    
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Flow Velocity (μm/s)")
    
    ax.legend()
    
    print("a = {:f}\nn = {:f}\nb = {:f}".format(fits1[0][0],
              fits1[0][1],fits2[0][0]))
    
    plt.tight_layout()
    
    return ax