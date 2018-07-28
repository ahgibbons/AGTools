# -*- coding: utf-8 -*-
"""
Created on Fri May 18 10:00:40 2018

@author: andrew
"""

import matplotlib.pyplot as plt
from matplotlib import lines, markers
from cycler import cycler

import numpy as np

monochrome_marker = (cycler('color',['k']) * cycler('linestyle',['-','--',':','-.']) * 
              cycler('marker',['^','x','+','s','.']))

monochrome = (cycler('color',['k',(0.4,0.4,0.4)]) * cycler('linestyle',['-','--','-.']))

default_cycle = plt.rcParams['axes.prop_cycle']


def plotCycleDec():
    def realDec(plot_func):
        def newplot(*args,ax=None,savepath=None,dpi=300,**kwargs):
            if ax is None:
                fig = plt.figure(figsize=(6,4),dpi=dpi)
                ax = fig.add_subplot(111)
        
            ax.set_prop_cycle(cycle)
        
            plot_func(*(ax, *args),**kwargs)
        
            if savepath:
                plt.savefig(savepath)
        
            return ax
        return newplot
    return realDec

def monochromeDecorator(plot_func):
    def newplot(*args,ax=None,savepath=None,dpi=300,**kwargs):
        if ax is None:
            fig = plt.figure(figsize=(6,4),dpi=dpi)
            ax = fig.add_subplot(111)
        
        ax.set_prop_cycle(monochrome)
        
        plot_func(*(ax, *args),**kwargs)
        
        if savepath:
            plt.savefig(savepath)
        
        return ax
    return newplot

def plotDecorator(plot_func):
    def newplot(*args,ax=None,savepath=None,dpi=300,**kwargs):
        if ax is None:
            fig = plt.figure(figsize=(6,4),dpi=dpi)
            ax = fig.add_subplot(111)
        
        plot_func(*(ax, *args),**kwargs)
        
        if savepath:
            plt.savefig(savepath)
        
        return ax
    return newplot


