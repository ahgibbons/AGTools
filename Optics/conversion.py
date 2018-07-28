# -*- coding: utf-8 -*-
"""
Created on Fri Jun 22 17:58:52 2018

@author: andrew
"""

import numpy as np

Gamma = 0.80
IntensityMax = 255

static private double Gamma = 0.80;
static private double IntensityMax = 255;


def waveLengthToRGB(wavelength):
    if (wavelength >= 380) and (wavelength<440):
        red   = -1*(wavelength - 440) / (440 - 380)
        green = 0.0
        blue  = 1.0
    elif (wavelength >= 440) and (wavelength<490):
        red   = 0.0
        green = (wavelength - 440) / (490 - 440)
        blue  = 1.0
    elif (wavelength >= 490) and (wavelength < 510):
        red   = 0.0
        green = 1.0
        blue = -1*(wavelength - 510) / (510 - 490)
    elif (wavelength >= 510) and (wavelength < 580):
        red   = (wavelength - 510) / (580 - 510)
        green = 1.0
        blue  = 0.0
    elif (wavelength >= 580) and (wavelength < 645):
        red = 1.0
        green = -1*(wavelength - 645) / (645 - 580)
        blue = 0.0
    elif (wavelength >= 645) and (wavelength < 781):
        red   = 1.0
        green = 0.0
        blue  = 0.0
    else:
        red   = 0
        green = 0
        blue  = 0

public static int[] waveLengthToRGB(double Wavelength){
    double factor;
    double Red,Green,Blue;

    if((Wavelength >= 380) && (Wavelength<440)){
        Red = -(Wavelength - 440) / (440 - 380);
        Green = 0.0;
        Blue = 1.0;
    }else if((Wavelength >= 440) && (Wavelength<490)){
        Red = 0.0;
        Green = (Wavelength - 440) / (490 - 440);
        Blue = 1.0;
    }else if((Wavelength >= 490) && (Wavelength<510)){
        Red = 0.0;
        Green = 1.0;
        Blue = -(Wavelength - 510) / (510 - 490);
    }else if((Wavelength >= 510) && (Wavelength<580)){
        Red = (Wavelength - 510) / (580 - 510);
        Green = 1.0;
        Blue = 0.0;
    }else if((Wavelength >= 580) && (Wavelength<645)){
        Red = 1.0;
        Green = -(Wavelength - 645) / (645 - 580);
        Blue = 0.0;
    }else if((Wavelength >= 645) && (Wavelength<781)){
        Red = 1.0;
        Green = 0.0;
        Blue = 0.0;
    }else{
        Red = 0.0;
        Green = 0.0;
        Blue = 0.0;
    };

    // Let the intensity fall off near the vision limits

    if((Wavelength >= 380) && (Wavelength<420)){
        factor = 0.3 + 0.7*(Wavelength - 380) / (420 - 380);
    }else if((Wavelength >= 420) && (Wavelength<701)){
        factor = 1.0;
    }else if((Wavelength >= 701) && (Wavelength<781)){
        factor = 0.3 + 0.7*(780 - Wavelength) / (780 - 700);
    }else{
        factor = 0.0;
    };


    int[] rgb = new int[3];

    // Don't want 0^x = 1 for x <> 0
    rgb[0] = Red==0.0 ? 0 : (int) Math.round(IntensityMax * Math.pow(Red * factor, Gamma));
    rgb[1] = Green==0.0 ? 0 : (int) Math.round(IntensityMax * Math.pow(Green * factor, Gamma));
    rgb[2] = Blue==0.0 ? 0 : (int) Math.round(IntensityMax * Math.pow(Blue * factor, Gamma));

    return rgb;
}