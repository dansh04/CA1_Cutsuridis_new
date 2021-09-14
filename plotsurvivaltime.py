#%%
import numpy as np
import matplotlib.pyplot as plt

pattern_cells = [1,2,7,11,21,28,35,38,39,43,46,49,56,57,59,62,78,81,88,90]
mgconcs = [0.0, 0.5, 1.0, 2.0, 4.0] 
pattern_avgtime = []
nonpattern_avgtime = []

for mgconc in mgconcs:
    data = np.loadtxt('pyresults/mgf_{}_cell_death.dat'.format(mgconc),skiprows=1)
    if len(data) != 0: 
        dead_times = data[:,0]
        dead_cells = data[:,1]        
    pcount = 0
    npcount = 0
    ptime = 0
    nptime = 0
    for i in range(len(data)):
        if dead_cells[i] in pattern_cells:
            pcount += 1
            ptime += dead_times[i]
        else:
            npcount += 1
            nptime += dead_times[i]
    for j in range(len(pattern_cells)-pcount): ptime += 1550       
    for k in range(100-len(pattern_cells)-npcount): nptime += 1550
    pattern_avgtime.append(ptime/len(pattern_cells))
    nonpattern_avgtime.append(nptime/(100-len(pattern_cells)))

#%%    
with open('pyresults/mgf_pattern_avgtime.dat', 'w') as f:
    f.write("mg\t time\n")
    for r in range(len(pattern_avgtime)):
        f.write("{}\t {}\n".format(mgconcs[r],pattern_avgtime[r]))

with open('pyresults/mgf_nonpattern_avgtime.dat', 'w') as f:
    f.write("mg\t time\n")
    for r in range(len(nonpattern_avgtime)):
        f.write("{}\t {}\n".format(mgconcs[r],nonpattern_avgtime[r]))

#%%        
x = np.arange(len(mgconcs))
plt.figure(figsize=(12,8))
plt.plot(x,nonpattern_avgtime,c='#2C599D',label='Non-Pattern',ls='-',lw=1.2,marker='s',markersize=6)
plt.plot(x,pattern_avgtime,c='#F98125',label='Pattern',ls='-',lw=1.2,marker='o',markersize=6)
plt.xticks(x,mgconcs,fontsize='large')
plt.yticks(fontsize='large')
plt.xlabel(xlabel='Magnesium Concentration (mM)',fontsize='x-large',fontweight='medium',loc='center',labelpad=14)
plt.ylabel(ylabel='Average Survival Time (ms)',fontsize='x-large',fontweight='medium',loc='center',labelpad=14)
plt.title(label='Average Survival Time vs. Mg Concentration',fontsize='xx-large',fontweight='roman',loc='center',pad=20)
plt.legend(loc='upper left',fontsize='large',fancybox=True,edgecolor='0.6')
plt.grid(linewidth=0.2)
plt.show()
