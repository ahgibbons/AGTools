# -*- coding: utf-8 -*-
"""
Created on Mon Jun 11 20:42:10 2018

@author: andrew
"""

import Sivaniah.Plot as Plot
import Sivaniah.Flow.flowspeed as flowspeed
import subprocess
import matplotlib.pyplot as plt
import os


test_file= r"G:\test\20180607 AG channel penta15C video 1 - 1 sec\20180607 AG channel penta15C video 1 - 1 sec_t001.AVI"
test_out =  r"G:\test\20180607 AG channel penta15C video 1 - 1 sec"
def extract_images(source_file,target_dir,scale=1):
    
    os.chdir(target_dir)
    
    subprocess.call(["ffmpeg","-i",test_file,"-vf","scale=iw/{:d}:ih/{:d}".format(scale,scale),
                    os.path.join(target_dir,"img_%04d.png")])
    

        
        