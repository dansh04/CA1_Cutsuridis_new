#%%
import numpy as np
import matplotlib.pyplot as plt

pattern_cells = [1,2,7,11,21,28,35,38,39,43,46,49,56,57,59,62,78,81,88,90]
mgconcs = [0.0, 0.5, 1.0, 2.0, 4.0] 
pattern_count = []
nonpattern_count = []

for mgconc in mgconcs:
    data = np.loadtxt('pyresults/mgf_{}_cell_death.dat'.format(mgconc),skiprows=1)
    if len(data) != 0: 
        dead_times = data[:,0]
        dead_cells = data[:,1]        
    pcount = 0
    npcount = 0
    for i in range(len(data)):
        if dead_cells[i] in pattern_cells: pcount += 1
        else: npcount += 1
    pattern_count.append(pcount)
    nonpattern_count.append(npcount)

#%%
with open('pyresults/mgf_pattern_count.dat', 'w') as f:
    f.write("mg\t count\n")
    for r in range(len(pattern_count)):
        f.write("{}\t {}\n".format(mgconcs[r],pattern_count[r]))

with open('pyresults/mgf_nonpattern_count.dat', 'w') as f:
    f.write("mg\t count\n")
    for r in range(len(nonpattern_count)):
        f.write("{}\t {}\n".format(mgconcs[r],nonpattern_count[r]))

#%%
width = 0.35
x = np.arange(len(mgconcs))
plt.figure(figsize=(12,8))
plt.bar(x-width/2,nonpattern_count,width,color='#2C599D',label='Nonpattern')
plt.bar(x+width/2,pattern_count,width,color='#F98125',label='Pattern',)
plt.xticks(x,mgconcs,fontsize='large')
plt.yticks(fontsize='large')
plt.xlabel(xlabel='Magnesium Concentration (mM)',fontsize='x-large',fontweight='medium',loc='center',labelpad=14)
plt.ylabel(ylabel='Number of Dead Cells',fontsize='x-large',fontweight='medium',loc='center',labelpad=14)
plt.title(label='Dead Cell Count vs. Mg Concentration',fontsize='xx-large',fontweight='roman',loc='center',pad=20)
plt.legend(loc='upper right',fontsize='large',fancybox=True,edgecolor='0.6')
plt.show()  
