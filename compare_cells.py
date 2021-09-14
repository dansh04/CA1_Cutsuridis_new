# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 14:40:04 2020

@author: mbezaire
"""
from neuron import h
import numpy as np
h.load_file("nrngui.hoc")
h.celsius = 34
# create one of each cell type
# olm_cell2, axoaxonic_cell17S, basket_cell17S, bistratified_cell13S, pyramidal_cell_14Vb.hoc
h.load_file("pyramidal_cell_14Vb.hoc")
h('objref hoccell')
h('hoccell = new PyramidalCell()')
print("Made a hoc cell")

from model import cellClasses

pycell = cellClasses.PyramidalCell()
print("Made a python cell")

# add a current clamp to each one
newclamp = h.IClamp(pycell.soma(0.5))
newclamp.delay = 10
newclamp.dur = 500
newclamp.amp = .1

h('objref oldclamp')
h('hoccell.soma oldclamp = new IClamp(0.5)')
h.oldclamp.delay = 10
h.oldclamp.dur = 500
h.oldclamp.amp = .1

h.tstop = 500
print("Set current injection")
#%%
# record membrane potential, run a simulation, and compare results
py_soma_v = h.Vector().record(pycell.soma(0.5)._ref_v)
hoc_soma_v = h.Vector().record(h.hoccell.soma(0.5)._ref_v)
print("Ready to record")

h.run()

timevec = np.arange(0,h.tstop+h.dt, h.dt)
print("Time to plot")

import matplotlib.pyplot as plt
plt.figure()
plt.plot(timevec,py_soma_v,'g-',timevec,hoc_soma_v,'r-')
plt.show()
#%%
cells = []
cells.append(pycell)
cells.append(h.hoccell)

celldict=[]
celldict.append(cells[0].soma.psection())
celldict.append(cells[1].soma.psection())
print("Things about the cell somata that seem different:")
for key in celldict[0]:
    if (celldict[0][key] != celldict[1][key]):
        print(key)
#%%
for key in celldict[0]["density_mechs"]:
    if (celldict[0]["density_mechs"][key] != celldict[1]["density_mechs"][key]):
        print(key)  
        print(celldict[0]["density_mechs"][key])
        print(celldict[1]["density_mechs"][key])
        print(" ")
#%%
secnamelist=[]
for sec in cells[0].all:
    tmp = sec.name()
    i = tmp.find(".")
    secnamelist.append(tmp[i+1:])

for sec in secnamelist: #,"radTprox","radTmed","radTdist","lm_thick2","lm_medium2","lm_thin2","lm_thick1","lm_medium1","lm_thin1","oriprox1",
#"oridist1","oriprox2","oridist2","axon"}:
    exec("secref0 = cells[0]."+sec)
    exec("secref1 = cells[1]."+sec)
    if (secref0.nseg != secref1.nseg):
        print(sec,secref0.nseg, secref1.nseg)
    if (secref0.L != secref1.L):
        print(sec,secref0.L, secref1.L)
    if (secref0.diam != secref1.diam):
        print(sec,secref0.diam, secrf1.diam)

