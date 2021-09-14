#%%
import numpy as np
import matplotlib.pyplot as plt

pattern = [1,2,7,11,21,28,35,38,39,43,46,49,56,57,59,62,78,81,88,90]
mgconcs = [0.0, 0.5, 1.0, 1.5, 2.0]
colors = ['#5f0000','#000a6a','#9fa000','#053100','#480064']

#%%
for gid in range(100):
    if gid not in pattern:
        plt.figure(figsize=(12,3*len(mgconcs)))
        for i, mgconc in enumerate(mgconcs):
            spks = np.loadtxt("pyresults/mgf_{}_cellv_{}.dat".format(mgconc,gid),skiprows=1)
            plt.subplot(len(mgconcs),1,i+1)
            plt.plot(spks[:,0],spks[:,1],color=colors[i],label='Mg = {} mM'.format(mgconc))
            plt.ylabel("Membrane Potential (mV)")
            plt.ylim([-120, 50])
            plt.legend(loc='upper right')
            if i == 0:
                if gid in pattern:
                    plt.title('Pattern Pyramidal Cell {}'.format(gid))
                else:
                    plt.title('Nonpattern Pyramidal Cell {}'.format(gid))
        plt.xlabel('Time (ms)')
        plt.show()
    
#%%
plt.figure(figsize=(14,3*len(mgconcs)))
for i, mgconc in enumerate(mgconcs):
    voltage = np.loadtxt('pyresults/mgf_{}_cellv_1.dat'.format(mgconc),skiprows=1)
    plt.subplot(len(mgconcs),1,i+1)
    plt.plot(voltage[:,0],voltage[:,1],c=colors[i],label='Mg = {} mM'.format(mgconc),ls='-',lw=1.6)
    plt.xticks(fontsize='medium')
    plt.yticks(fontsize='medium')
    plt.ylabel(ylabel='Membrane Potential (mV)',fontsize='large',fontweight='medium',loc='center',labelpad=14)
    plt.xlim([-50,1600])
    plt.ylim([-115,45])
    plt.legend(loc='upper right',fontsize='medium',fancybox=True,edgecolor='0.6')
    if i==0:
        plt.title('Pattern Pyramidal Cell',fontsize='xx-large',fontweight='roman',loc='center',pad=20)
plt.xlabel(xlabel='Time (ms)',fontsize='large',fontweight='medium',loc='center',labelpad=14)
plt.show()

#%%
plt.figure(figsize=(14,3*len(mgconcs)))
for i, mgconc in enumerate(mgconcs):
    voltage = np.loadtxt('pyresults/mgf_{}_cellv_74.dat'.format(mgconc),skiprows=1)
    plt.subplot(len(mgconcs),1,i+1)
    plt.plot(voltage[:,0],voltage[:,1],c=colors[i],label='Mg = {} mM'.format(mgconc),ls='-',lw=1.6)
    plt.xticks(fontsize='medium')
    plt.yticks(fontsize='medium')
    plt.ylabel(ylabel='Membrane Potential (mV)',fontsize='large',fontweight='medium',loc='center',labelpad=14)
    plt.xlim([-50,1600])
    plt.ylim([-115,45])
    plt.legend(loc='upper right',fontsize='medium',fancybox=True,edgecolor='0.6')
    if i==0:
        plt.title('Non-Pattern Pyramidal Cell',fontsize='xx-large',fontweight='roman',loc='center',pad=20)
plt.xlabel(xlabel='Time (ms)',fontsize='large',fontweight='medium',loc='center',labelpad=14)
plt.show()
