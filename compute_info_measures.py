# -*- coding: utf-8 -*-
"""
Created on Thu Jul  5 18:24:45 2018


@author: Sam
"""

import numpy as np

def compute_info_measures(pdf):
    #pdf obtained from numpy histogramdd() function for a 2D pdf
    Hx1 = np.zeros(1)
    Hx2 = np.zeros(1)
    Ix1x2 = np.zeros(1)
    #marginal values (sum along rows and columns)
    m_i = np.sum(pdf,axis=1)
    m_j = np.sum(pdf,axis=0)
    #entropy values
    H1_vect = -1 * m_i * np.log2(m_i)
    #H1_vect[np.isnan(H1_vect)]=0
    Hx1 = np.sum(H1_vect)
    H2_vect = -1 * m_j * np.log2(m_j)
    H2_vect[np.isnan(H2_vect)]=0
    Hx2 = np.sum(H2_vect)
    #local information measures
    pi_pj = np.outer(m_i, m_j)
    f_ij = pdf/pi_pj   
    f_ij[f_ij == 0]=10**-10
    logf = np.where(pdf != 0, np.log2(f_ij), 0)
    plogf = pdf*logf
    #the sum of all plogf values = mutual informaton
    Ix0 = np.sum(plogf[0])
    Ix1 = np.sum(plogf[1])
    Ix1x2 = np.sum(plogf)
    #calculate specific information
    SI = (Ix0/Ix1)*(m_i[1]/m_i[0])
    infodict = {'Hx1':Hx1, 'Hx2':Hx2, 'Ix0':Ix0, 'Ix1':Ix1, 'Ix1x2':Ix1x2,
                'f':f_ij,'logf':logf,'plogf':plogf, 'SI':SI, 'm1':m_i[1], 'm0':m_i[0]}
    return infodict
