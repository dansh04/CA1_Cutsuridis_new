# Dummy cell containing a BurstStim object
# BPG 10-12-08
from neuron import h

class StimCell():  
    """ Stim cell with stim attribute that references a RegnStim point process """
    def __init__(self, gid=-1):
        self.gid = gid
        self.is_art =1
        self.ncstim = []
        #self.stim = h.BurstStim2()
        self.stim = h.RegnStim()
       	self.stim.number = 10000
       	self.stim.start = 0
       	self.stim.interval = 25
       	self.stim.noise = 0
           
    def __repr__(self):
        return "Stim cell {}: ISI is {} ms".format(self.gid, self.stim.interval)

    def connect2target(self,target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.ncstim.append(h.NetCon(self.stim, target))
        self.ncstim[-1].delay = delay # ms
        self.ncstim[-1].weight[0] = weight # NetCon weight is a vector    
        return self.ncstim[-1]