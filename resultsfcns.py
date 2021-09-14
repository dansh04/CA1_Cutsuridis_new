# Plot spikes from list of time vs cell number
# and determine quality of recall
# Summary diagram for Hippocampus paper
# BPG 23-1-09
import numpy as np
import math
import matplotlib.pyplot as plt

def plot_results(fstem,scaleDown=1,NUMCYCLES=8):
    scaleDown=1
    
    NCELL = 235-(1-scaleDown)*230  # number of cells (neurons); CA3, EC, SEP Pyr can be scaled down (230)
    NPCELL = 100*scaleDown # number of PC (output) cells
    SPATT = 20*scaleDown   # number of active cells per pattern
    
    
    if scaleDown<1:
        FPATT = r'Weights\pattsN100S20P5Scaled.dat' # TODO: Replace with your full path to the file
    else:
        FPATT = r'Weights\pattsN100S20P5.dat' # TODO: Replace with your full path to the file
        
    NPATT = 1   # number of patterns
    CPATT = 1  # index of cue pattern
    
    RTIME = 50+(250*NUMCYCLES)    # run time (msecs)
    STIME = 200
    ETIME = 2050
    
    
    patts = np.loadtxt(FPATT)   # load stored patterns
    cue = patts[:,CPATT]   # extract cue pattern
    
    FSPIKE = r'pyresults/{}_spt.dat'.format(fstem)   # spikes file
    sp = np.loadtxt(FSPIKE,skiprows=1)  # load spike times
    st = sp[:,0]       # extract times
    cell = sp[:,1]     # extract corresponding cell indices
    # extract PC spiking
    stp = [x for i, x in enumerate(st) if cell[i]<NPCELL ]
    #stp = st(cell < NPCELL)
    cellp = [x for i, x in enumerate(cell) if cell[i]<NPCELL ]
    
    #cellp = cell(cell < NPCELL)
    
    # Analyse spiking over time and compare with cue
    DT = 1 # sliding time
    #TW = 5    # width of sliding time window
    TW = 10    # width of sliding time window
    
    ti = range(0,RTIME-TW,DT)
    NW = len(ti)   # number of time windows
    nc = np.zeros((NW,1))
    ha = np.zeros((NW,1))
    co = np.zeros((NW,1))
    an = np.zeros((NW,1))
    mc = cue.mean() # mean cue activity
    
    for i in range(NW):
        # TODO: Pythonize the matlab line below:
        # rp = cellp(stp>=ti[i] & stp<ti[i]+TW) # active cells in sliding window
        # My first attempt at the line above resulted in an IndexError: list index out of range:
        rp = [cellp[int(x)] for x in stp if (x>ti[i] and x<(ti[i]+TW))]  # active cells in sliding window
    
        nc[i] = len(rp)    # number of active cells in window
    
        p = np.zeros((NPCELL))
        
        for x in rp:
            p[int(x)+1] = 1 #     p[rp+1,1] = 1  # recalled pattern
            
        ha[i] = (sum(p == cue)/NPCELL)  # hamming distance
        mp = p.mean()   # mean pattern activity
        # correlation (normalised dot product)
        if mp == 0:
            co[i] = 0
        else:
            # TODO check the next line to convert from MATLAB to Python syntax
            co[i] = np.dot(p,cue)/math.sqrt(sum(p)*sum(cue))
    
        # angle (Graham & Willshaw 97)
        if mp == 0 or mp == 1:
            an[i] = 0
        else:
            # TODO check the next line to convert from MATLAB to Python syntax
            an[i] = sum((p-mp)*(cue-mc))/math.sqrt(sum((p-mp)**2)*sum((cue-mc)**2))
    
    
    # Generate figure
    plt.figure()
    ms=8
    lw=2
    
    plt.subplot(4,1,1)
    plt.plot(sp[:,0], sp[:,1], 'k.', marker='o')   # raster plot of Sep, EC & CA3 spiking
    plt.title('(a) Input spikes')
    plt.ylabel('Cell no.')
    plt.xlim([STIME, ETIME])
    plt.ylim([NPCELL+4, NPCELL+4+130])
    plt.subplot(4,1,2)
    plt.plot(sp[:,0], sp[:,1], 'k.', marker='o')   # raster plot of PC spiking
    plt.title('(b) Pyramidal cell spikes')
    plt.ylabel('Cell no.')
    plt.xlim([STIME, ETIME])
    plt.ylim([0, NPCELL-1])
    
    plt.subplot(4,1,3)
    #hold on
    plt.plot(ti, nc, 'k-', LineWidth=2) # spike counts
    plt.title('(c) PC spike count')
    plt.ylabel('Spike count')
    plt.xlim([STIME, ETIME])
    plt.ylim([0, NPCELL])
    plt.subplot(4,1,4)
    #hold on
    plt.plot(ti, co, 'k-', LineWidth=2) # recall quality
    plt.title('(d) Recall quality')
    plt.ylabel('Quality')
    plt.xlabel('Time (msecs)')
    plt.xlim([STIME, ETIME])
    plt.ylim([0, 1.02])
    
    plt.savefig("Images/{}.png".format(fstem))
    plt.show()
    
    print(co[co>0].mean())


plot_results("myname")    