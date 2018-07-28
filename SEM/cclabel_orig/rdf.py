# -*- coding: utf-8 -*-
"""
Created on Wed May 16 12:59:21 2018

@author: andrew
"""

## Radial Distribution Analysis

import cclabel
from PIL import Image
from PIL import ImageFilter
from PIL import ImageOps
from math import sqrt

area = 1280*1024

def imgPillarSpacing(imgpath,imgdims=(1280,1024),scale=4.673,threshold=150):
    img = Image.open(imgpath)
    w,h = imgdims
    img_crop = img.crop((0,0,w,h)).convert()
    print("Image cropped...")
    blurred_img = img_crop.filter(ImageFilter.GaussianBlur(radius=2)).convert()
    print("Gaussian blur applied...")
    th_img = blurred_img.point(lambda p: p>threshold and 255).convert()
    th_img = ImageOps.invert(th_img).convert('1')
    print("Threshold applied to image...")
    labels = cclabel.ccRun(th_img)
    print("Connected regions calculated...")
    value_set = set(labels.values())
    print("{:d} regions...".format(len(value_set)))
    #coms = [labelCentre(labels,i) for i in set(labels.values())]
    N = len(value_set)
    sigma = N/area
    d = 1/sqrt(sigma)
    print("Spacing is {:d}.".format(d))
    print("Finished.\n")
    
    return d*scale
    
    

def labelCentre(ccdata,index):
    x_t = 0
    y_t = 0
    n = 0
    
    for k in ccdata.keys():
        v = ccdata[k]
        if v==index:
            x_t += k[0]
            y_t += k[1]
            n += 1
    x = x_t / n
    y = y_t / n
    return (x,y)
    