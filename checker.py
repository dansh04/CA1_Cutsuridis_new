# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 10:27:32 2020

@author: mbezaire
"""


import numpy as np

vals = np.loadtxt("Results/HAM_P5R1_AAC.dat")
pyvals = np.loadtxt("pyresults/HAM_P5R1_AAC.dat")

import matplotlib.pyplot as plt

plt.figure()
plt.plot(pyvals, "y-", vals, "k--")
plt.xlim([0, 5000])
plt.show()