# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 15:47:43 2021

@author: maria
"""

from neuron import h
import numpy as np


h.load_file("stdrun.hoc")
h.load_file("nrngui.hoc") # load_file


class BallAndStick:
    def __init__(self, gid):
        self._gid = gid
        self.soma = h.Section(name='soma', cell=self)
        self.dend = h.Section(name='dend', cell=self)
        self.dend.connect(self.soma)
        self.soma.L = self.soma.diam = 12.6157
        # self.syn = h.NMDA2(self.soma(0.5))
        # self.syn.gmax = 0.00065
        # self.syn.mg = 1
        # self.syn.onset = 10
        # self.syn.tauF = 9
        # self.syn.tauS = 90

        self.syn = h.NMDA(self.soma(0.5))
        self.syn.gNMDAmax = 1
        self.syn.tcon = 2.3    
        self.syn.tcoff = 100
       
        # self.syn = h.MyExp2Syn(self.soma(0.5)) 
        # self.syn.tau1 = 0.5
        # self.syn.tau2 = 3  
        
stim = h.NetStim()
stim.number = 1
stim.start = 5
stim.interval = 0
stim.noise = 0

my_cell = BallAndStick(0)

cells = [ my_cell, stim]

nc = h.NetCon(stim, my_cell.syn)
nc.delay = 3
nc.weight[0] = .005

h.tstop = 500
syn_current = h.Vector().record(my_cell.syn._ref_i)
soma_v = h.Vector().record(my_cell.soma(0.5)._ref_v)

tvec = h.Vector()
idvec = h.Vector()

nc0 = h.NetCon(stim, None)
nc0.record(tvec, idvec, 0)
nc0.delay = 3
nc0.weight[0] = 0.004

nc1 = h.NetCon(my_cell.soma(0.5)._ref_v, None, sec=my_cell.soma)
nc1.record(tvec, idvec, 1)
nc1.threshold = -10
nc1.delay = 3
nc1.weight[0] = 0.004      
        
h.run()

import matplotlib.pyplot as plt

t = np.arange(0,h.tstop+h.dt,h.dt)
i = []

plt.plot(t,soma_v)
plt.ylabel('Membrane Potential (mV)')
plt.figure()
plt.plot(t,syn_current)
plt.ylabel('Synaptic Current (nA)')
