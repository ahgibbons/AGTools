# -*- coding: utf-8 -*-
"""
Created on Mon Sep  3 12:20:17 2018

@author: andrew
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


class Ray:
    def __init__(self, origin, direction):
        self.origin = np.array(origin)
        self.direction = np.array(direction)
        
        
class Sphere:
    def __init__(self,origin,radius):
        self.origin = np.array(origin)
        self.radius = radius
        

def normal_vector_of_incidence(ray, sphere):
    intersection_points = intersect_line_sphere(ray, sphere)
    intersection_point = intersection_points[1]
    
    normal = intersection_point - sphere.origin
    normal = normal / np.sqrt(np.dot(normal,normal))
    
    return normal
    
def snells_law(incident,index1,index2):
    sin_i = np.sin(incident)
    sin_r = index1/index2 * sin_i
    return np.arcsin(sin_r)
    
def intersect_line_sphere(ray,sphere):
     c = sphere.origin
     r_sq = sphere.radius**2
     o = ray.origin
     l = ray.direction
     
     l_sq = np.dot(l,l)
     o_sq = np.dot(o,o)
     c_sq = np.dot(c,c)
     
     determinant = np.sqrt((np.dot(c-o,l))**2 - l_sq*(o_sq + c_sq - 2*np.dot(c,o) - r_sq))
     
     lambda1 = (np.dot(c-o,l) + determinant) / l_sq
     
     lambda2 = (np.dot(c-o,l) - determinant) / l_sq
     
     p1 = o + lambda1*l
     p2 = o + lambda2*l
     
     return (p1, p2)
 
ray = Ray(np.array([0,0,0]),np.array([-1,0,0]))
ray2 = Ray([0,0.5,0],[-1,0,0])
sphere = Sphere(np.array([0,0,0]),1)