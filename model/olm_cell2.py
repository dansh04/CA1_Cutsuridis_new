# Data from Saraga et al. (2003) paper
# changed morphology and some channel densities (BPG 12-1-09)
#   OLM_Cell
# basic_shape,
# pre_list, connect2target
# soma, dend1, dend2, axon
# all
from neuron import h

class OLMCell():
    """ OLM Cell definition """
    def __init__(self, gid = -1):
        self.x = 0; self.y = 0; self.z = 0
        self.gid = gid
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        # pre_list = new List()
        self.addSynapses() # synapses
        self.is_art = 0
        self.nc = []
        
    def __repr__(self):
        return "OLM Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend1 = h.Section(name='dend1', cell=self)
        self.dend2 = h.Section(name='dend2', cell=self)
        self.axon = h.Section(name='axon', cell=self)

    def build_topology(self):
        self.dend1.connect(self.soma(1))
        self.dend2.connect(self.soma(0))
        self.axon.connect(self.soma(1))
        # TODO basic_shape()

    def define_geometry(self):
        self.soma.pt3dclear()
        self.dend1.pt3dclear()
        self.dend2.pt3dclear()
        self.axon.pt3dclear()
        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 3, sec=self.dend1)
        h.pt3dadd(90, 0, 0, 3, sec=self.dend1)
        h.pt3dadd(0, 0, 0, 3, sec=self.dend2)
        h.pt3dadd(-74, 0, 0, 3, sec=self.dend2)
        h.pt3dadd(15, 0, 0, 1.5, sec=self.axon)
        h.pt3dadd(15, 150, 0, 1.5, sec=self.axon)
        self.soma.L = 20
        self.soma.diam = 10
        self.dend1.L = 250
        self.dend1.diam = 3
        self.dend2.L = 250
        self.dend2.diam = 3
        self.axon.L = 150
        self.axon.diam = 1.5
        
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)


    def define_biophysics(self):
        self.Rm = 20000 # 1/5e-05    
        
        for sec in self.all:
            self.Ra = 150
            self.cm = 1.3

        self.soma.insert("IA")
        self.soma.insert("Ih")
        self.soma.insert("Ksoma")
        self.soma.insert("Nasoma")
        for seg in self.soma:
            seg.gkAbar_IA = 0.0165
            seg.gkhbar_Ih = 0.0005 # 0.001385
            seg.gksoma_Ksoma = 0.0319
            seg.gnasoma_Nasoma = 0.0107
            seg.gl_Nasoma = 1/self.Rm
            seg.el_Nasoma = -70

        self.dend1.insert("IA")
        # self.dend1.insert("Ih")
        self.dend1.insert("Kdend")
        self.dend1.insert("Nadend")
        for seg in self.dend1:
            seg.gkAbar_IA = 0.004 # 0.013
            # seg.gkhbar_Ih = 0.0005 # 0.001385
            seg.gkdend_Kdend = 2*0.023
            seg.gnadend_Nadend = 2*0.0117
            seg.gl_Nadend = 1/self.Rm
            seg.el_Nadend = -70

        self.dend2.insert("IA")
        # self.dend2.insert("Ih")
        self.dend2.insert("Kdend")
        self.dend2.insert("Nadend")
        for seg in self.dend2:
            seg.gkAbar_IA = 0.004 # 0.013
            # seg.gkhbar_Ih = 0.0005 # 0.001385
            seg.gkdend_Kdend = 2*0.023
            seg.gnadend_Nadend = 2*0.0117
            seg.gl_Nadend = 1/self.Rm
            seg.el_Nadend = -70


        self.axon.insert("Kaxon")
        self.axon.insert("Naaxon")
        for seg in self.axon:
            seg.gkaxon_Kaxon = 0.05104
            seg.gnaaxon_Naaxon = 0.01712
            seg.gl_Naaxon = 1/self.Rm
            seg.el_Naaxon = -70


    def connect2target(self,target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]


    def addSynapses(self):
        self.pre_list = []
        
        # E0
        syn_ = h.Exp2Syn(self.dend2(0.5))
        self.pre_list.append(syn_)    # AMPA        Pyramidal
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E1
        syn_ = h.Exp2Syn(self.dend1(0.5))
        self.pre_list.append(syn_)    # AMPA        Pyramidal
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I2
        syn_ = h.Exp2Syn(self.dend1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I3
        syn_ = h.Exp2Syn(self.dend1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

