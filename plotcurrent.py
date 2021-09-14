#%% 
import numpy as np
import matplotlib.pyplot as plt
    
mgconcs = [round(0.2*i,1) for i in list(range(11))]
avgcurrents = [-0.02, -0.015, -0.01, -0.005, -0]
nonpattern_avgcur = []
pattern_avgcur = []

for mgconc in mgconcs: 
    nonpattern_spks = np.loadtxt("pyresults/mgf_{}_celli_0.dat".format(mgconc),skiprows=1)
    nonpattern_cot = nonpattern_spks[:,1]
    nonpattern_avgcur.append(np.mean(nonpattern_cot))
    pattern_spks = np.loadtxt("pyresults/mgf_{}_celli_1.dat".format(mgconc),skiprows=1)
    pattern_cot = pattern_spks[:,1]
    pattern_avgcur.append(np.mean(pattern_cot))

#%%
x = np.arange(len(mgconcs))
plt.figure(figsize=(12,8))
plt.plot(x,nonpattern_avgcur,c='#2C599D',label='Non-Pattern',ls='-',lw=1.2,marker='s',markersize=6)
plt.plot(x,pattern_avgcur,c='#F98125',label='Pattern',ls='-',lw=1.2,marker='o',markersize=6)
plt.xticks(x,mgconcs,fontsize='large')
plt.yticks(fontsize='large')
plt.xlabel(xlabel='Magnesium Concentration (mM)',fontsize='x-large',fontweight='medium',loc='center',labelpad=14)
plt.ylabel(ylabel='Average Current (nA)',fontsize='x-large',fontweight='medium',loc='center',labelpad=14)
plt.title(label='Average Current vs. Mg Concentration',fontsize='xx-large',fontweight='roman',loc='center',pad=20)
plt.legend(loc='upper left',fontsize='large',fancybox=True,edgecolor='0.6')
plt.grid(linewidth=0.2)
plt.show()
