# -*- coding: utf-8 -*-
"""
Compute Icit.py
Created on Thu Jul  5 18:09:48 2018

@author: Sam
"""

import numpy as np
from compute_info_measures import compute_info_measures as cim

#compute I critical values to compare to I values to test significance
def compute_icrit(pptlist_new, sflist_new):  
    nTests = 100
    xbin = 2
    ybin = 5
    Itot_vector = []
    for j in range(0, nTests):
        ppt_shuffled = list(pptlist_new)
        np.random.shuffle(ppt_shuffled)                        
        Tuple = np.array((ppt_shuffled, sflist_new)).T
        pdf, edges = np.histogramdd(Tuple, (xbin, ybin))
        pdf = pdf/np.sum(pdf)
        info = cim(pdf)
        Itot_vector.append(info['Ix1x2'])
        I_mean=np.average(Itot_vector)
        I_std=np.std(Itot_vector) 
        Icrit = I_mean + 3*I_std
        Icritdict = {'Icrit':Icrit}
        return Icritdict
