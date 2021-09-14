# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 10:30:39 2020

@author: mbezaire
"""

runs2include=["par","sim1"]

performance_list=[]
death_list=[]
electrostim_list=[]

for run in runs2include:
    with open('pyresults/' + run+'_performance.txt', 'r') as f:  # Python 3: open(..., 'wb')
        content = f.readlines()
    
        perf = float(content[0])
        electrostim = float(content[1])
        percentDeath = float(content[2])

        performance_list.append(perf)
        death_list.append(electrostim)
        electrostim_list.append(percentDeath)
        
import matplotlib.pyplot as plt

plt.figure()
plt.plot(death_list,performance_list,'ro')
plt.xlabel('% Cell Death')
plt.ylabel('Memory Recall Performance (Scale of 0 to 1)')
plt.show()
plt.savefig('Images/memory_recall_v_death.png')

