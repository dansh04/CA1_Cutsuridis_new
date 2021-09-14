# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 11:00:21 2020

@author: mbezaire
"""


import numpy as np
import matplotlib.pyplot as plt
rez={'npvsoma':0,
 'npvsr': 0,
 'npvslm': 0,
 'pvsoma':0,
 'pvsr': 0,
 'pvslm': 0,
 'BC': 1,
 'AAC': 1,
 'BSC': 1,
 'OLM': 1}

for key in rez:
    if rez[key]==1:
        newway = np.loadtxt(r"C:\Users\mbezaire\Documents\GitHub\CA1_Cutsuridis\pyresults\par_v"+key+"_0.dat")
        #oldway = np.loadtxt(r"C:\Users\mbezaire\Documents\GitHub\CA1_Cutsuridis\pyresults\check_v"+key+"_0.dat")
    else:
        newway = np.loadtxt(r"C:\Users\mbezaire\Documents\GitHub\CA1_Cutsuridis\pyresults\par_"+key+"_0.dat")
        #oldway = np.loadtxt(r"C:\Users\mbezaire\Documents\GitHub\CA1_Cutsuridis\pyresults\check_"+key+"_0.dat")
    oldway = np.loadtxt(r"C:\Users\mbezaire\Desktop\123815-master\Results\HAM_P5R1_"+key+".dat")

    plt.figure()
    plt.plot(newway,'k-', oldway, 'r--')
    plt.title(key)
    plt.show()