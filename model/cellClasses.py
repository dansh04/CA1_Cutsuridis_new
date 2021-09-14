# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 21:53:18 2020

@author: mbezaire
"""
keepoldtypo = 0

from neuron import h

class modelcell():
    def __init__(self):
        self.x = 0; self.y = 0; self.z = 0

        self.gid = -1
        self.core_i=-1
        self.coretype_i=-1
        
        self.is_art = 0
        self.nc = []
        self.pre_list = []

class PyramidalCell(modelcell):
    """ Pyramidal Cell definition """
    def __init__(self, gid = -1, mgconc = 1):
        super().__init__()
        self.gid = gid
        self.mgconc = mgconc
        self.list_syns = []
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        self.stimobj = h.IClamp(self.soma(0.5))
        self.stimobj.delay = 1e9
        self.stimobj.dur = 1
        self.stimobj.amp = 0  
        self.isdead = 0
        # self.addSynapses() # synapses
        
    def __repr__(self):
        return "Pyramidal Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radTprox = h.Section(name='radTprox', cell=self)
        self.radTmed = h.Section(name='radTmed', cell=self)
        self.radTdist = h.Section(name='radTdist', cell=self)
        self.axon = h.Section(name='axon', cell=self)
        self.lm_thick2 = h.Section(name='lm_thick2', cell=self)
        self.lm_thick1 = h.Section(name='lm_thick1', cell=self)
        self.oriprox1 = h.Section(name='oriprox1', cell=self)
        self.oridist1 = h.Section(name='oridist1', cell=self)
        self.oriprox2 = h.Section(name='oriprox2', cell=self)
        self.oridist2 = h.Section(name='oridist2', cell=self)
        self.lm_medium2 = h.Section(name='lm_medium2', cell=self)
        self.lm_thin2 = h.Section(name='lm_thin2', cell=self)
        self.lm_medium1 = h.Section(name='lm_medium1', cell=self)
        self.lm_thin1 = h.Section(name='lm_thin1', cell=self)

    def build_topology(self):
        self.radTprox.connect(self.soma(1))
        self.radTmed.connect(self.radTprox(1))
        self.radTdist.connect(self.radTmed(1))
        self.axon.connect(self.soma(1))
        self.lm_thick2.connect(self.radTdist(1))
        self.lm_medium2.connect(self.lm_thick2(1))
        self.lm_thin2.connect(self.lm_medium2(1))

        self.lm_thick1.connect(self.radTdist(1))
        self.lm_medium1.connect(self.lm_thick1(1))
        self.lm_thin1.connect(self.lm_medium1(1))

        self.oriprox1.connect(self.soma(0))
        self.oridist1.connect(self.oriprox1(1))

        self.oriprox2.connect(self.soma(1))
        self.oridist2.connect(self.oriprox2(1))

    def define_geometry(self):
        for sec in self.all:
            sec.pt3dclear()

        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        
        h.pt3dadd(15, 0, 0, 1, sec=self.radTprox)
        h.pt3dadd(15, 30, 0, 1, sec=self.radTprox)
       
        h.pt3dadd(15, 30, 0, 1, sec=self.radTmed)
        h.pt3dadd(15, 60, 0, 1, sec=self.radTmed)
        
        h.pt3dadd(15, 60, 0, 1, sec=self.radTdist)
        h.pt3dadd(15, 90, 0, 1, sec=self.radTdist)

        h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick2)
        h.pt3dadd(45, 105, 0, 1, sec=self.lm_thick2)

        h.pt3dadd(45, 105, 0, 1, sec=self.lm_medium2)
        h.pt3dadd(75, 120, 0, 1, sec=self.lm_medium2)

        h.pt3dadd(75, 120, 0, 1, sec=self.lm_thin2)
        h.pt3dadd(105, 135, 0, 1, sec=self.lm_thin2)

        h.pt3dadd(15, 90, 0, 1, sec=self.lm_thick1)
        h.pt3dadd(-14, 105, 0, 1, sec=self.lm_thick1)

        h.pt3dadd(-14, 105, 0, 1, sec=self.lm_medium1)
        h.pt3dadd(-44, 120, 0, 1, sec=self.lm_medium1)

        h.pt3dadd(-44, 120, 0, 1, sec=self.lm_thin1)
        h.pt3dadd(-89, 135, 0, 1, sec=self.lm_thin1)
         
        h.pt3dadd(15, 0, 0, 1, sec=self.axon)
        h.pt3dadd(15, -149, 0, 1, sec=self.axon)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriprox1)
        h.pt3dadd(-44, -29, 0, 1, sec=self.oriprox1)

        h.pt3dadd(-44, -29, 0, 1, sec=self.oridist1)
        h.pt3dadd(-74, -59, 0, 1, sec=self.oridist1)

        h.pt3dadd(15, 0, 0, 1, sec=self.oriprox2)
        h.pt3dadd(60, -29, 0, 1, sec=self.oriprox2)

        h.pt3dadd(60, -29, 0, 1, sec=self.oridist2)
        h.pt3dadd(105, -59, 0, 1, sec=self.oridist2)
    
        self.soma.L = 10
        self.soma.diam = 10
        self.radTprox.L = 100
        self.radTprox.diam = 4
        self.radTmed.L = 100
        self.radTmed.diam = 3
        self.radTdist.L = 200
        self.radTdist.diam = 2

        self.lm_medium2.L = 100
        self.lm_medium2.diam = 1.5
        self.lm_thin2.L = 50
        self.lm_thin2.diam = 1
        self.lm_thick2.L = 100
        self.lm_thick2.diam = 2
        
        self.lm_thick1.L = 100
        self.lm_thick1.diam = 2
        self.lm_medium1.L = 100
        self.lm_medium1.diam = 1.5
        self.lm_thin1.L = 50
        self.lm_thin1.diam = 1

        self.oriprox1.L = 100
        self.oriprox1.diam = 2
        self.oridist1.L = 200
        self.oridist1.diam = 1.5

        self.oriprox2.L = 100
        self.oriprox2.diam = 2
        self.oridist2.L = 200
        self.oridist2.diam = 1.5

        self.axon.L = 150
        self.axon.diam = 1
 
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 # TODO This equation is the same but sometimes produces smaller numbers of nseg than the hoc one

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)

    def define_biophysics(self):
        gka_soma = 0.0075
        gh_soma = 0.00005
        Rm = 20000    # 28000 Ohm.cm^2 (Migliore value)
        
        self.soma.insert("hha2") # HH mechanism with low threshold for Na spikes (-57 mV)
        self.soma.insert("pas") # leak conductance
        self.soma.insert("h") # h current according to Migliore et al. 2004 
        self.soma.insert("kap") # proximal A current
        self.soma.insert("km") # m-type potassium current
        self.soma.insert("cal") # HVA Ca++-L type current
        self.soma.insert("cat") # LVA Ca++-T type current
        self.soma.insert("somacar") # HVAm Ca++-R type current
        self.soma.insert("kca") # K(Ca) sAHP potassium type current
        self.soma.insert("mykca") # medium AHP K++ current (BPG)
        self.soma.insert("cad") # calcium pump/buffering mechanism
        
        for seg in self.soma:
            seg.gnabar_hha2 = 0.007
            seg.gkbar_hha2 = 0.007/5
            seg.gl_hha2 = 0
            seg.el_hha2 = -70
            seg.g_pas =  1/Rm

            seg.ghdbar_h = gh_soma
            seg.vhalfl_h = -73
            
            seg.gkabar_kap = gka_soma # 0.0075
            seg.gbar_km = 0.06
            seg.gcalbar_cal = 0.0014/2
            seg.gcatbar_cat = 0.0001/2
            seg.gcabar_somacar = 0.0003
            seg.gbar_kca = 5*0.0001
            seg.gkbar_mykca = 0.09075
  
  # //        insert hNa              // h current according to Poirazi 2003
  # //        gbar_h  = 0.000043      // anything above 0.000043 gives hyperpolarizing oscillations
  # //        gbar_h  = 1.872e-5        
  # //        K_h     = 8.8
  # //        vhalf_h = -82
  
        self.radTprox.insert("h") # h current according to Migliore et al. 2004 
        self.radTprox.insert("car") # h current according to Migliore et al. 2004 
        self.radTprox.insert("calH") # h current according to Migliore et al. 2004 
        self.radTprox.insert("cat") # h current according to Migliore et al. 2004 
        self.radTprox.insert("cad") # calcium pump/buffering mechanism
        self.radTprox.insert("kca") # slow AHP K+ current
        self.radTprox.insert("mykca") # medium AHP K++ current (BPG)
        self.radTprox.insert("km") # m-type K current
        self.radTprox.insert("kap") # Inserting A-current
        self.radTprox.insert("kad") # Inserting A-current
        self.radTprox.insert("hha_old") # HH mechanism with high threshold for Na spikes (-50 mV)
        self.radTprox.insert("pas") # passive
        
        for seg in self.radTprox:
            seg.ghdbar_h = 2*gh_soma # 0.000005    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001
            seg.gbar_kca = 5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 2*gka_soma # 0.0075
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70
        
  # //        insert hNa             // h current according to Poirazi 2003
  # //        gbar_h  = 0.000043     // anything above 0.000043 gives hyperpolarizing oscillations
  # //        gbar_h  = 1.872e-5        
  # //        K_h     = 8.8
  # //        vhalf_h = -82

        self.radTmed.insert("h") # h current according to Migliore et al. 2004
        self.radTmed.insert("car") # HVAm Ca++-R type current
        self.radTmed.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.radTmed.insert("cat") # HVA T-type Ca2+ channel 
        self.radTmed.insert("cad") # calcium pump/buffering mechanism
        self.radTmed.insert("kca") # slow AHP K+ current
        self.radTmed.insert("mykca") # medium AHP K++ current (BPG)
        self.radTmed.insert("km") # m-type K current
        self.radTmed.insert("kap") # Inserting A-current
        self.radTmed.insert("kad") 
        self.radTmed.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.radTmed.insert("pas") # leak conductance
        
        for seg in self.radTmed:
            seg.ghdbar_h = 4*gh_soma # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 10*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001 # 0.0001
            seg.gbar_kca = 5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0
            seg.gkabar_kad = 4*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70
            
  # //        insert hNa              // h current according to Poirazi 2003
  # //        gbar_h  = 0.000043        
  # //        gbar_h  = 1.872e-5        
  # //        K_h     = 8.8
  # //        vhalf_h = -82

        self.radTdist.insert("h") # h current according to Migliore et al. 2004
        self.radTdist.insert("car") # HVAm Ca++-R type current
        self.radTdist.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.radTdist.insert("cat") # HVA T-type Ca2+ channel 
        self.radTdist.insert("cad") # calcium pump/buffering mechanism
        self.radTdist.insert("kca") # slow AHP K+ current
        self.radTdist.insert("mykca") # medium AHP K++ current (BPG)
        self.radTdist.insert("km") # m-type K current
        self.radTdist.insert("kap") # Inserting A-current
        self.radTdist.insert("kad") 
        self.radTdist.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.radTdist.insert("pas") # leak conductance
        
        for seg in self.radTdist:
            seg.ghdbar_h = 7*gh_soma # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 10*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001 # 0.0001
            seg.gbar_kca = 0.5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 0.25*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0
            seg.gkabar_kad = 6*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70

        self.lm_thick2.insert("kad") 
        self.lm_thick2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thick2.insert("pas") # leak conductance
        
        for seg in self.lm_thick2:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000

        self.lm_medium2.insert("kad") 
        self.lm_medium2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_medium2.insert("pas") # leak conductance
        
        for seg in self.lm_medium2:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000

        self.lm_thin2.insert("kad") 
        self.lm_thin2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thin2.insert("pas") # leak conductance
        
        for seg in self.lm_thin2:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000

        self.lm_thick1.insert("kad") 
        self.lm_thick1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thick1.insert("pas") # leak conductance
        
        for seg in self.lm_thick1:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000

        self.lm_medium1.insert("kad") 
        self.lm_medium1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_medium1.insert("pas") # leak conductance
        
        for seg in self.lm_medium1:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000

        self.lm_thin1.insert("kad") 
        self.lm_thin1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.lm_thin1.insert("pas") # leak conductance
        
        for seg in self.lm_thin1:
            seg.gkabar_kad = 6.5*gka_soma
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old  = 0.007/8.065
            seg.el_hha_old = -70
            seg.g_pas = 1/200000

        self.oriprox1.insert("h") # h current according to Migliore et al. 2004
        self.oriprox1.insert("car") # HVAm Ca++-R type current
        self.oriprox1.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oriprox1.insert("cat") # HVA T-type Ca2+ channel 
        self.oriprox1.insert("cad") # calcium pump/buffering mechanism
        self.oriprox1.insert("kca") # slow AHP K+ current
        self.oriprox1.insert("mykca") # medium AHP K++ current (BPG)
        self.oriprox1.insert("km") # m-type K current
        self.oriprox1.insert("kap") # Inserting A-current
        self.oriprox1.insert("kad") 
        self.oriprox1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oriprox1.insert("pas") # leak conductance
        
        for seg in self.oriprox1:
            seg.ghdbar_h = gh_soma # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001 # 0.0001
            seg.gbar_kca = 5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = gka_soma
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70
            
        self.oriprox2.insert("h") # h current according to Migliore et al. 2004
        self.oriprox2.insert("car") # HVAm Ca++-R type current
        self.oriprox2.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oriprox2.insert("cat") # HVA T-type Ca2+ channel 
        self.oriprox2.insert("cad") # calcium pump/buffering mechanism
        self.oriprox2.insert("kca") # slow AHP K+ current
        self.oriprox2.insert("mykca") # medium AHP K++ current (BPG)
        self.oriprox2.insert("km") # m-type K current
        self.oriprox2.insert("kap") # Inserting A-current
        self.oriprox2.insert("kad") 
        self.oriprox2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oriprox2.insert("pas") # leak conductance
        
        for seg in self.oriprox2:
            seg.ghdbar_h = gh_soma # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001 # 0.0001
            seg.gbar_kca = 5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0.0075
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70

        self.oridist1.insert("h") # h current according to Migliore et al. 2004
        self.oridist1.insert("car") # HVAm Ca++-R type current
        self.oridist1.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oridist1.insert("cat") # HVA T-type Ca2+ channel 
        self.oridist1.insert("cad") # calcium pump/buffering mechanism
        self.oridist1.insert("kca") # slow AHP K+ current
        self.oridist1.insert("mykca") # medium AHP K++ current (BPG)
        self.oridist1.insert("km") # m-type K current
        self.oridist1.insert("kap") # Inserting A-current
        self.oridist1.insert("kad") 
        self.oridist1.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oridist1.insert("pas") # leak conductance
        
        for seg in self.oridist1:
            seg.ghdbar_h = 2*gh_soma # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001 # 0.0001
            seg.gbar_kca = 5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = gka_soma
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70

        self.oridist2.insert("h") # h current according to Migliore et al. 2004
        self.oridist2.insert("car") # HVAm Ca++-R type current
        self.oridist2.insert("calH") # HVA L-type Ca2+ channel used in distal dendrites to account for distally restricted initiation of Ca2+ spikes
        self.oridist2.insert("cat") # HVA T-type Ca2+ channel 
        self.oridist2.insert("cad") # calcium pump/buffering mechanism
        self.oridist2.insert("kca") # slow AHP K+ current
        self.oridist2.insert("mykca") # medium AHP K++ current (BPG)
        self.oridist2.insert("km") # m-type K current
        self.oridist2.insert("kap") # Inserting A-current
        self.oridist2.insert("kad") 
        self.oridist2.insert("hha_old") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.oridist2.insert("pas") # leak conductance
        
        for seg in self.oridist2:
            seg.ghdbar_h = 2*gh_soma # 0.000005                    
            seg.vhalfl_h = -81
            seg.gcabar_car = 0.1*0.0003
            seg.gcalbar_calH = 0.1*0.00031635 # 4.6*0.00031635 varies from .1*0.00031635 to 4.6*0.00031635 as distance increases
            seg.gcatbar_cat = 0.0001 # 0.0001
            seg.gbar_kca = 5*0.0001 # varies depending on distance from 0.5*0.0001 to 5*0.0001
            seg.gkbar_mykca = 2*0.0165
            seg.gbar_km = 0.06 # varies with distance (see Poirazzi et al. 2003 cell-setup.hoc file)
            seg.gkabar_kap = 0.0075
            seg.gkabar_kad = 0
            seg.gnabar_hha_old = 0.007
            seg.gkbar_hha_old = 0.007/8.065
            seg.el_hha_old = -70

        self.axon.insert("km") 
        self.axon.insert("hha2") # // HH mechanism with high threshold for Na spikes (-50 mV)
        self.axon.insert("pas") # leak conductance
        
        for seg in self.axon:
            seg.gbar_km = 0.5*0.06
            seg.gnabar_hha2 = 0.1
            seg.gkbar_hha2 = 0.1/5
            seg.gl_hha2 = 0
            seg.el_hha2 = -70
            seg.g_pas = 1/Rm

        for sec in self.all:
            # self.cm = Not setting cm
            sec.Ra = 150 # 31.3 +/- 10.9
            sec.cm = 1
            sec.ena = 50
            sec.e_pas = -70
            sec.g_pas = 1/Rm # crucial parameter for backpropagating action potential spiking of PCs
            sec.ek = -80

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]

    def addSynapses(self):
        self.pre_list = []

        # E0
        syn_ = h.MyExp2Syn(self.lm_thick1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 0.5
        syn_.e = 0
        
        # E1
        syn_ = h.MyExp2Syn(self.lm_thick2(0.5))
        self.pre_list.append(syn_)    # AMPA        EC 
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E2
        syn_ = h.MyExp2Syn(self.radTmed(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E3
        syn_ = h.NMDA2(self.radTmed(0.5))
        self.pre_list.append(syn_)    # NMDA        CA3 Shaffer collateral
        syn_.gmax = 0.00065
        syn_.mg = self.mgconc
        # syn_.tcon = 2.3    
        # syn_.tcoff = 100
        # syn_.gNMDAmax = 1    # use connection weight to determine max cond

        # E4
        syn_ = h.MyExp2Syn(self.radTprox(0.5))
        self.pre_list.append(syn_)    # AMPA        PC Recurrent collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I5
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A    basket cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I6
        syn_ = h.MyExp2Syn(self.axon(0.1))
        self.pre_list.append(syn_)    # GABA-A    AA cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I7
        syn_ = h.MyExp2Syn(self.lm_thick1(0.5))
        self.pre_list.append(syn_)    # GABA-A    OLM cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I8
        syn_ = h.MyExp2Syn(self.lm_thick2(0.5))
        self.pre_list.append(syn_)    # GABA-A    OLM cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I9
        syn_ = h.MyExp2Syn(self.lm_thick1(0.5))
        self.pre_list.append(syn_)    # GABA-B    OLM cell
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # I10
        syn_ = h.MyExp2Syn(self.lm_thick2(0.5))
        self.pre_list.append(syn_)    # GABA-B    OLM Cell
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # I11
        syn_ = h.MyExp2Syn(self.radTmed(0.8))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I12
        syn_ = h.MyExp2Syn(self.radTmed(0.7))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I13
        syn_ = h.MyExp2Syn(self.radTmed(0.6))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
                
        # I14
        syn_ = h.MyExp2Syn(self.radTmed(0.4))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
                
        # I15
        syn_ = h.MyExp2Syn(self.radTmed(0.3))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75       
                
        # I16
        syn_ = h.MyExp2Syn(self.radTmed(0.2))
        self.pre_list.append(syn_)    # GABA-A    Bistratified
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75

        # I17
        syn_ = h.MyExp2Syn(self.radTmed(0.8))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
                
        # I18
        syn_ = h.MyExp2Syn(self.radTmed(0.7))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
                
        # I19
        syn_ = h.MyExp2Syn(self.radTmed(0.6))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75        
                
        # I20
        syn_ = h.MyExp2Syn(self.radTmed(0.4))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

        # I21
        syn_ = h.MyExp2Syn(self.radTmed(0.3))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75       
                
        # I22
        syn_ = h.MyExp2Syn(self.radTmed(0.2))
        self.pre_list.append(syn_)    # GABA-B    Bistratified
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # E23
        syn_ = h.STDPE2(self.radTmed(0.5))
        self.pre_list.append(syn_)    # AMPA modifiable	CA3 Schaffer collaterals
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
class OLMCell(modelcell):
    """ OLM Cell definition """
    def __init__(self, gid = -1):
        super().__init__()
        self.gid = gid
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        self.addSynapses() # synapses

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
        
        for sec in self.all:
            sec.Ra = 150
            sec.cm = 1.3

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
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

class BasketCell(modelcell):
    """ Basket Cell definition """
    def __init__(self, gid = -1):
        super().__init__()
        self.gid = gid
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        self.addSynapses() # synapses
        
    def __repr__(self):
        return "Basket Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radT2 = h.Section(name='radT2', cell=self)
        self.radM2 = h.Section(name='radM2', cell=self)
        self.radt2 = h.Section(name='radt2', cell=self)
        self.radT1 = h.Section(name='radT1', cell=self)
        self.radM1 = h.Section(name='radM1', cell=self)
        self.radt1 = h.Section(name='radt1', cell=self)
        self.oriT1 = h.Section(name='oriT1', cell=self)
        self.oriM1 = h.Section(name='oriM1', cell=self)
        self.orit1 = h.Section(name='orit1', cell=self)
        self.oriT2 = h.Section(name='oriT2', cell=self)
        self.oriM2 = h.Section(name='oriM2', cell=self)
        self.orit2 = h.Section(name='orit2', cell=self)
        self.lmM2 = h.Section(name='lmM2', cell=self)
        self.lmt2 = h.Section(name='lmt2', cell=self)
        self.lmM1 = h.Section(name='lmM1', cell=self)
        self.lmt1 = h.Section(name='lmt1', cell=self)

    def build_topology(self):
        self.radT2.connect(self.soma(1))
        self.radM2.connect(self.radT2(1))
        self.radt2.connect(self.radM2(1))
        self.radT1.connect(self.soma(0))
        self.radM1.connect(self.radT1(1))
        self.radt1.connect(self.radM1(1))

        self.oriT1.connect(self.soma(0))
        self.oriM1.connect(self.oriT1(1))
        self.orit1.connect(self.oriM1(1))

        self.oriT2.connect(self.soma(1))
        self.oriM2.connect(self.oriT2(1))
        self.orit2.connect(self.oriM2(1))

        self.lmM2.connect(self.radt2(1))
        self.lmt2.connect(self.lmM2(1))

        self.lmM1.connect(self.radt1(1))
        self.lmt1.connect(self.lmM1(1))

    def define_geometry(self):
        for sec in self.all:
            sec.pt3dclear()

        h.pt3dadd(-44, 45, 0, 1, sec=self.lmM1)
        h.pt3dadd(-59, 60, 0, 1, sec=self.lmM1)

        h.pt3dadd(-59, 60, 0, 1, sec=self.lmt1)
        h.pt3dadd(-89, 90, 0, 1, sec=self.lmt1)

        h.pt3dadd(105, 90, 0, 1, sec=self.lmt2)
        h.pt3dadd(120, 105, 0, 1, sec=self.lmt2)

        h.pt3dadd(90, 75, 0, 1, sec=self.lmM2)
        h.pt3dadd(105, 90, 0, 1, sec=self.lmM2)

        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        
        h.pt3dadd(15, 0, 0, 1, sec=self.radT2)
        h.pt3dadd(45, 30, 0, 1, sec=self.radT2)
        
        h.pt3dadd(45, 30, 0, 1, sec=self.radM2)
        h.pt3dadd(75, 60, 0, 1, sec=self.radM2)
        
        h.pt3dadd(75, 60, 0, 1, sec=self.radt2)
        h.pt3dadd(90, 75, 0, 1, sec=self.radt2)
        
        h.pt3dadd(0, 0, 0, 1, sec=self.radT1)
        h.pt3dadd(-14, 15, 0, 1, sec=self.radT1)

        h.pt3dadd(-14, 15, 0, 1, sec=self.radM1)
        h.pt3dadd(-29, 30, 0, 1, sec=self.radM1)

        h.pt3dadd(-29, 30, 0, 1, sec=self.radt1)
        h.pt3dadd(-44, 45, 0, 1, sec=self.radt1)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriT1)
        h.pt3dadd(-29, -29, 0, 1, sec=self.oriT1)

        h.pt3dadd(-29, -29, 0, 1, sec=self.oriM1)
        h.pt3dadd(-59, -59, 0, 1, sec=self.oriM1)

        h.pt3dadd(-59, -59, 0, 1, sec=self.orit1)
        h.pt3dadd(-89, -89, 0, 1, sec=self.orit1)

        h.pt3dadd(15, 0, 0, 1, sec=self.oriT2)
        h.pt3dadd(45, -29, 0, 1, sec=self.oriT2)

        h.pt3dadd(45, -29, 0, 1, sec=self.oriM2)
        h.pt3dadd(75, -59, 0, 1, sec=self.oriM2)

        h.pt3dadd(75, -59, 0, 1, sec=self.orit2)
        h.pt3dadd(105, -89, 0, 1, sec=self.orit2)
       
        self.soma.L = 20
        self.soma.diam = 10
        self.radT2.L = 100
        self.radT2.diam = 4
        self.radM2.L = 100
        self.radM2.diam = 3
        self.radt2.L = 200
        self.radt2.diam = 2

        self.lmM2.L = 100
        self.lmM2.diam = 1.5
        self.lmt2.L = 100
        self.lmt2.diam = 1

        self.radT1.L = 100
        self.radT1.diam = 4
        self.radM1.L = 100
        self.radM1.diam = 3
        self.radt1.L = 200
        self.radt1.diam = 2

        self.oriT1.L = 100
        self.oriT1.diam = 2
        self.oriM1.L = 100
        self.oriM1.diam = 1.5
        self.orit1.L = 100
        self.orit1.diam = 1

        self.oriT2.L = 100
        self.oriT2.diam = 2
        self.oriM2.L = 100
        self.oriM2.diam = 1.5
        self.orit2.L = 100
        self.orit2.diam = 1

        self.lmM1.L = 100
        self.lmM1.diam = 1.5
        self.lmt1.L = 100
        self.lmt1.diam = 1
 
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)

    def define_biophysics(self):
        gna = 0.2
        
        self.soma.insert("ichan2")
        for seg in self.soma:
            seg.gnatbar_ichan2 = gna  		# 0.12 //original 0.030 to .055 
            seg.gkfbar_ichan2 = 0.013  		# original 0.015
            seg.gl_ichan2 = 0.00018
        self.soma.cm=1.4
        
        self.lmt1.insert("ichan2")
        for seg in self.lmt1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmt1.cm=1.4
        
        self.lmt2.insert("ichan2")
        for seg in self.lmt2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmt2.cm=1.4
        
        self.lmM1.insert("ichan2")
        for seg in self.lmM1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmM1.cm=1.4
        
        self.lmM2.insert("ichan2")
        for seg in self.lmM2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmM2.cm=1.4
        
        self.radt1.insert("ichan2")
        for seg in self.radt1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radt1.cm=1.4
        
        self.radt2.insert("ichan2")
        for seg in self.radt2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radt2.cm=1.4
        
        self.radM1.insert("ichan2")
        for seg in self.radM1:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radM1.cm=1.4
        
        self.radM2.insert("ichan2")
        for seg in self.radM2:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radM2.cm=1.4
        
        self.radT1.insert("ichan2")
        for seg in self.radT1:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radT1.cm=1.4
        
        self.radT2.insert("ichan2")
        for seg in self.radT2:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radT2.cm=1.4
        
        self.oriT1.insert("ichan2")
        for seg in self.oriT1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriT1.cm=1.4
        
        self.oriT2.insert("ichan2")
        for seg in self.oriT2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriT2.cm=1.4
        
        self.oriM1.insert("ichan2")
        for seg in self.oriM1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriM1.cm=1.4
        
        self.oriM2.insert("ichan2")
        for seg in self.oriM2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriM2.cm=1.4
        
        self.orit1.insert("ichan2")
        for seg in self.orit1:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
        self.orit1.cm=1.4
        
        self.orit2.insert("ichan2")
        for seg in self.orit2:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
        self.orit2.cm=1.4
        
        for sec in self.all:            
            sec.insert("ccanl")
            for seg in sec:
                seg.catau_ccanl = 10		# Time constant for decay of intracellular Ca2+
                seg.caiinf_ccanl = 5.e-6		# Steady-state intracellular Ca2+ concentration

            sec.insert("borgka")
            for seg in sec:
                seg.gkabar_borgka = 0.00015		# A-type K+ conductance
            
            sec.insert("nca") # N-type Ca2+ conductance
            for seg in sec:
                seg.gncabar_nca = 0.0008   		# check to modify- original 0.004
           
            sec.insert("lca") 
            for seg in sec:
                seg.glcabar_lca = 0.005		# L-type Ca2+ conductance
           
            sec.insert("gskch") 
            for seg in sec:
                seg.gskbar_gskch = 0.000002		# Ca2+-dependent K (SK) conductance
           
            sec.insert("mykca") 
            if (keepoldtypo==0):
                # commenting out gbar to match their typo
                for seg in sec:
                    seg.gkbar_mykca = 0.0002			# Ca2+ and Voltage-dependent K+ (BK) conductance
		 					# make catau slower70e-3 	cao=2 cai=50.e-6
            # self.cm = Not setting cm
            sec.Ra = 100			# 31.3 +/- 10.9
            sec.enat = 55
            sec.ekf = -90
            sec.ek = -90
            sec.elca = 130
            sec.esk = -90
            sec.el_ichan2 = -60			#-60.06
        self.cao_ccanl = 2

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]

    def addSynapses(self):
        self.pre_list = []

        # E0
        syn_ = h.MyExp2Syn(self.lmM1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E1
        syn_ = h.MyExp2Syn(self.lmM2(0.5))
        self.pre_list.append(syn_)    # AMPA        EC 
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E2
        syn_ = h.MyExp2Syn(self.radM1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E3
        syn_ = h.MyExp2Syn(self.radM2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E4
        syn_ = h.MyExp2Syn(self.radT1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E5
        syn_ = h.MyExp2Syn(self.radT2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E6
        syn_ = h.MyExp2Syn(self.oriT1(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E7
        syn_ = h.MyExp2Syn(self.oriT2(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I8
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A	Neighboring basket cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I9
        syn_ = h.MyExp2Syn(self.soma(0.6))
        self.pre_list.append(syn_)    # GABA-A	Bistratified cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I10
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I11
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I12
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

        # I13
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

class AACell(modelcell):
    """ Axo-axonic Cell definition """
    def __init__(self, gid = -1):
        super().__init__()
        self.gid = gid
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        self.addSynapses() # synapses
        
    def __repr__(self):
        return "Axo-axonic Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radT2 = h.Section(name='radT2', cell=self)
        self.radM2 = h.Section(name='radM2', cell=self)
        self.radt2 = h.Section(name='radt2', cell=self)
        self.radT1 = h.Section(name='radT1', cell=self)
        self.radM1 = h.Section(name='radM1', cell=self)
        self.radt1 = h.Section(name='radt1', cell=self)
        self.oriT1 = h.Section(name='oriT1', cell=self)
        self.oriM1 = h.Section(name='oriM1', cell=self)
        self.orit1 = h.Section(name='orit1', cell=self)
        self.oriT2 = h.Section(name='oriT2', cell=self)
        self.oriM2 = h.Section(name='oriM2', cell=self)
        self.orit2 = h.Section(name='orit2', cell=self)
        self.lmM2 = h.Section(name='lmM2', cell=self)
        self.lmt2 = h.Section(name='lmt2', cell=self)
        self.lmM1 = h.Section(name='lmM1', cell=self)
        self.lmt1 = h.Section(name='lmt1', cell=self)

    def build_topology(self):
        self.radT2.connect(self.soma(1))
        self.radM2.connect(self.radT2(1))
        self.radt2.connect(self.radM2(1))
        self.radT1.connect(self.soma(0))
        self.radM1.connect(self.radT1(1))
        self.radt1.connect(self.radM1(1))

        self.oriT1.connect(self.soma(0))
        self.oriM1.connect(self.oriT1(1))
        self.orit1.connect(self.oriM1(1))

        self.oriT2.connect(self.soma(1))
        self.oriM2.connect(self.oriT2(1))
        self.orit2.connect(self.oriM2(1))

        self.lmM2.connect(self.radt2(1))
        self.lmt2.connect(self.lmM2(1))

        self.lmM1.connect(self.radt1(1))
        self.lmt1.connect(self.lmM1(1))

    def define_geometry(self):
        for sec in self.all:
            sec.pt3dclear()

        h.pt3dadd(-44, 45, 0, 1, sec=self.lmM1)
        h.pt3dadd(-59, 60, 0, 1, sec=self.lmM1)

        h.pt3dadd(-59, 60, 0, 1, sec=self.lmt1)
        h.pt3dadd(-89, 90, 0, 1, sec=self.lmt1)

        h.pt3dadd(105, 90, 0, 1, sec=self.lmt2)
        h.pt3dadd(120, 105, 0, 1, sec=self.lmt2)

        h.pt3dadd(90, 75, 0, 1, sec=self.lmM2)
        h.pt3dadd(105, 90, 0, 1, sec=self.lmM2)

        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        
        h.pt3dadd(15, 0, 0, 1, sec=self.radT2)
        h.pt3dadd(45, 30, 0, 1, sec=self.radT2)
        
        h.pt3dadd(45, 30, 0, 1, sec=self.radM2)
        h.pt3dadd(75, 60, 0, 1, sec=self.radM2)
        
        h.pt3dadd(75, 60, 0, 1, sec=self.radt2)
        h.pt3dadd(90, 75, 0, 1, sec=self.radt2)
        
        h.pt3dadd(0, 0, 0, 1, sec=self.radT1)
        h.pt3dadd(-14, 15, 0, 1, sec=self.radT1)

        h.pt3dadd(-14, 15, 0, 1, sec=self.radM1)
        h.pt3dadd(-29, 30, 0, 1, sec=self.radM1)

        h.pt3dadd(-29, 30, 0, 1, sec=self.radt1)
        h.pt3dadd(-44, 45, 0, 1, sec=self.radt1)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriT1)
        h.pt3dadd(-29, -29, 0, 1, sec=self.oriT1)

        h.pt3dadd(-29, -29, 0, 1, sec=self.oriM1)
        h.pt3dadd(-59, -59, 0, 1, sec=self.oriM1)

        h.pt3dadd(-59, -59, 0, 1, sec=self.orit1)
        h.pt3dadd(-89, -89, 0, 1, sec=self.orit1)

        h.pt3dadd(15, 0, 0, 1, sec=self.oriT2)
        h.pt3dadd(45, -29, 0, 1, sec=self.oriT2)

        h.pt3dadd(45, -29, 0, 1, sec=self.oriM2)
        h.pt3dadd(75, -59, 0, 1, sec=self.oriM2)

        h.pt3dadd(75, -59, 0, 1, sec=self.orit2)
        h.pt3dadd(105, -89, 0, 1, sec=self.orit2)
        
        self.soma.L = 20
        self.soma.diam = 10
        self.radT2.L = 100
        self.radT2.diam = 4
        self.radM2.L = 100
        self.radM2.diam = 3
        self.radt2.L = 200
        self.radt2.diam = 2

        self.lmM2.L = 100
        self.lmM2.diam = 1.5
        self.lmt2.L = 100
        self.lmt2.diam = 1

        self.radT1.L = 100
        self.radT1.diam = 4
        self.radM1.L = 100
        self.radM1.diam = 3
        self.radt1.L = 200
        self.radt1.diam = 2

        self.oriT1.L = 100
        self.oriT1.diam = 2
        self.oriM1.L = 100
        self.oriM1.diam = 1.5
        self.orit1.L = 100
        self.orit1.diam = 1

        self.oriT2.L = 100
        self.oriT2.diam = 2
        self.oriM2.L = 100
        self.oriM2.diam = 1.5
        self.orit2.L = 100
        self.orit2.diam = 1

        self.lmM1.L = 100
        self.lmM1.diam = 1.5
        self.lmt1.L = 100
        self.lmt1.diam = 1
 
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)

    def define_biophysics(self):
        gna = 0.15
        
        self.soma.insert("ichan2")
        for seg in self.soma:
            seg.gnatbar_ichan2 = gna  		# 0.12 //original 0.030 to .055 
            seg.gkfbar_ichan2 = 0.013  		# original 0.015
            seg.gl_ichan2 = 0.00018
        self.soma.cm=1.4
        
        self.lmt1.insert("ichan2")
        for seg in self.lmt1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmt1.cm=1.4
        
        self.lmt2.insert("ichan2")
        for seg in self.lmt2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmt2.cm=1.4
        
        self.lmM1.insert("ichan2")
        for seg in self.lmM1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmM1.cm=1.4
        
        self.lmM2.insert("ichan2")
        for seg in self.lmM2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.lmM2.cm=1.4
        
        self.radt1.insert("ichan2")
        for seg in self.radt1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radt1.cm=1.4
        
        self.radt2.insert("ichan2")
        for seg in self.radt2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radt2.cm=1.4
        
        self.radM1.insert("ichan2")
        for seg in self.radM1:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radM1.cm=1.4
        
        self.radM2.insert("ichan2")
        for seg in self.radM2:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radM2.cm=1.4
        
        self.radT1.insert("ichan2")
        for seg in self.radT1:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radT1.cm=1.4
        
        self.radT2.insert("ichan2")
        for seg in self.radT2:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radT2.cm=1.4
        
        self.oriT1.insert("ichan2")
        for seg in self.oriT1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriT1.cm=1.4
        
        self.oriT2.insert("ichan2")
        for seg in self.oriT2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriT2.cm=1.4

        self.oriM1.insert("ichan2")
        for seg in self.oriM1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriM1.cm=1.4
        
        self.oriM2.insert("ichan2")
        for seg in self.oriM2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriM2.cm=1.4
        
        self.orit1.insert("ichan2")
        for seg in self.orit1:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
        self.orit1.cm=1.4
        
        self.orit2.insert("ichan2")
        for seg in self.orit2:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
        self.orit2.cm=1.4
        
        for sec in self.all:
            # self.cm = Not setting cm
            
            sec.insert("ccanl")
            for seg in sec:
                seg.catau_ccanl = 10		# Time constant for decay of intracellular Ca2+
                seg.caiinf_ccanl = 5.e-6		# Steady-state intracellular Ca2+ concentration
            
            sec.insert("borgka")
            for seg in sec:
                seg.gkabar_borgka = 0.00015		# A-type K+ conductance
            
            sec.insert("nca") # N-type Ca2+ conductance
            for seg in sec:
                seg.gncabar_nca = 0.0008   		# check to modify- original 0.004
           
            sec.insert("lca") 
            for seg in sec:
                seg.glcabar_lca = 0.005		# L-type Ca2+ conductance
           
            sec.insert("gskch") 
            for seg in sec:
                seg.gskbar_gskch = 0.000002		# Ca2+-dependent K (SK) conductance
           
            sec.insert("mykca") 
            if (keepoldtypo==0):
                # commenting out gbar to match their typo
                for seg in sec:
                    seg.gkbar_mykca = 0.0002			# Ca2+ and Voltage-dependent K+ (BK) conductance
		 					# make catau slower70e-3 	cao=2 cai=50.e-6
            sec.Ra = 100			# 31.3 +/- 10.9
            sec.enat = 55
            sec.ekf = -90
            sec.ek = -90
            sec.elca = 130
            sec.esk = -90
            sec.el_ichan2 = -60			#-60.0
        self.cao_ccanl = 2

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]

    def addSynapses(self):
        self.pre_list = []

        # E0
        syn_ = h.MyExp2Syn(self.lmM1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E1
        syn_ = h.MyExp2Syn(self.lmM2(0.5))
        self.pre_list.append(syn_)    # AMPA        EC (not used)
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # E2
        syn_ = h.MyExp2Syn(self.radM1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E3
        syn_ = h.MyExp2Syn(self.radM2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E4
        syn_ = h.MyExp2Syn(self.radT1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E5
        syn_ = h.MyExp2Syn(self.radT2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E6
        syn_ = h.MyExp2Syn(self.oriT1(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E7
        syn_ = h.MyExp2Syn(self.oriT2(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I8
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A	Neighboring axo-axonic cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I9
        syn_ = h.MyExp2Syn(self.soma(0.6))
        self.pre_list.append(syn_)    # GABA-A	Bistratified cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I10
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I11
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I12
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # I13
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
class BistratifiedCell(modelcell):
    """ Bistratified Cell definition """
    def __init__(self, gid = -1):
        super().__init__()
        self.gid = gid
        self.create_sections() 
        self.build_topology()
        self.build_subsets() # subsets()
        self.define_geometry() # geom()
        self.define_biophysics() # biophys()
        self.addSynapses() # synapses
        
    def __repr__(self):
        return "Bistratified Cell {}".format(self.gid)

    def create_sections(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radT2 = h.Section(name='radT2', cell=self)
        self.radM2 = h.Section(name='radM2', cell=self)
        self.radt2 = h.Section(name='radt2', cell=self)
        self.radT1 = h.Section(name='radT1', cell=self)
        self.radM1 = h.Section(name='radM1', cell=self)
        self.radt1 = h.Section(name='radt1', cell=self)
        self.oriT1 = h.Section(name='oriT1', cell=self)
        self.oriM1 = h.Section(name='oriM1', cell=self)
        self.orit1 = h.Section(name='orit1', cell=self)
        self.oriT2 = h.Section(name='oriT2', cell=self)
        self.oriM2 = h.Section(name='oriM2', cell=self)
        self.orit2 = h.Section(name='orit2', cell=self)

    def build_topology(self):
        self.radT2.connect(self.soma(1))
        self.radM2.connect(self.radT2(1))
        self.radt2.connect(self.radM2(1))
        self.radT1.connect(self.soma(0))
        self.radM1.connect(self.radT1(1))
        self.radt1.connect(self.radM1(1))

        self.oriT1.connect(self.soma(0))
        self.oriM1.connect(self.oriT1(1))
        self.orit1.connect(self.oriM1(1))

        self.oriT2.connect(self.soma(1))
        self.oriM2.connect(self.oriT2(1))
        self.orit2.connect(self.oriM2(1))

    def define_geometry(self):
        for sec in self.all:
            sec.pt3dclear()

        h.pt3dadd(0, 0, 0, 10, sec=self.soma)
        h.pt3dadd(15, 0, 0, 10, sec=self.soma)
        
        h.pt3dadd(15, 0, 0, 1, sec=self.radT2)
        h.pt3dadd(45, 30, 0, 1, sec=self.radT2)
        
        h.pt3dadd(45, 30, 0, 1, sec=self.radM2)
        h.pt3dadd(75, 60, 0, 1, sec=self.radM2)
        
        h.pt3dadd(75, 60, 0, 1, sec=self.radt2)
        h.pt3dadd(90, 75, 0, 1, sec=self.radt2)
        
        h.pt3dadd(0, 0, 0, 1, sec=self.radT1)
        h.pt3dadd(-14, 15, 0, 1, sec=self.radT1)

        h.pt3dadd(-14, 15, 0, 1, sec=self.radM1)
        h.pt3dadd(-29, 30, 0, 1, sec=self.radM1)

        h.pt3dadd(-29, 30, 0, 1, sec=self.radt1)
        h.pt3dadd(-44, 45, 0, 1, sec=self.radt1)

        h.pt3dadd(0, 0, 0, 1, sec=self.oriT1)
        h.pt3dadd(-29, -29, 0, 1, sec=self.oriT1)

        h.pt3dadd(-29, -29, 0, 1, sec=self.oriM1)
        h.pt3dadd(-59, -59, 0, 1, sec=self.oriM1)

        h.pt3dadd(-59, -59, 0, 1, sec=self.orit1)
        h.pt3dadd(-89, -89, 0, 1, sec=self.orit1)

        h.pt3dadd(15, 0, 0, 1, sec=self.oriT2)
        h.pt3dadd(45, -29, 0, 1, sec=self.oriT2)

        h.pt3dadd(45, -29, 0, 1, sec=self.oriM2)
        h.pt3dadd(75, -59, 0, 1, sec=self.oriM2)

        h.pt3dadd(75, -59, 0, 1, sec=self.orit2)
        h.pt3dadd(105, -89, 0, 1, sec=self.orit2)
        
        self.soma.L = 20
        self.soma.diam = 10
        self.radT2.L = 100
        self.radT2.diam = 4
        self.radM2.L = 100
        self.radM2.diam = 3
        self.radt2.L = 200
        self.radt2.diam = 2

        self.radT1.L = 100
        self.radT1.diam = 4
        self.radM1.L = 100
        self.radM1.diam = 3
        self.radt1.L = 200
        self.radt1.diam = 2

        self.oriT1.L = 100
        self.oriT1.diam = 2
        self.oriM1.L = 100
        self.oriM1.diam = 1.5
        self.orit1.L = 100
        self.orit1.diam = 1

        self.oriT2.L = 100
        self.oriT2.diam = 2
        self.oriM2.L = 100
        self.oriM2.diam = 1.5
        self.orit2.L = 100
        self.orit2.diam = 1
        
        for sec in self.all:
            h("lf = lambda_f(100)")
            sec.nseg = int((sec.L/(0.1*h.lf)+.9)/2)*2 + 1 

    def build_subsets(self):
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)

    def define_biophysics(self):
        gna = 0.3
        
        self.soma.insert("ichan2")
        for seg in self.soma:
            seg.gnatbar_ichan2 = gna  		# 0.12 //original 0.030 to .055 
            seg.gkfbar_ichan2 = 0.013  		# original 0.015
            seg.gl_ichan2 = 0.00018
        self.soma.cm=1.4
        
        self.radt1.insert("ichan2")
        for seg in self.radt1:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radt1.cm=1.4
        
        self.radt2.insert("ichan2")
        for seg in self.radt2:
            seg.gnatbar_ichan2 = gna  		# 0.4  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radt2.cm=1.4
        
        self.radM1.insert("ichan2")
        for seg in self.radM1:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radM1.cm=1.4
        
        self.radM2.insert("ichan2")
        for seg in self.radM2:
            seg.gnatbar_ichan2 = gna  		# 0.3  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radM2.cm=1.4
        
        self.radT1.insert("ichan2")
        for seg in self.radT1:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radT1.cm=1.4
        
        self.radT2.insert("ichan2")
        for seg in self.radT2:
            seg.gnatbar_ichan2 = gna  		# 0.2  	//original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.radT2.cm=1.4
        
        self.oriT1.insert("ichan2")
        for seg in self.oriT1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriT1.cm=1.4
        
        self.oriT2.insert("ichan2")
        for seg in self.oriT2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriT2.cm=1.4
        
        self.oriM1.insert("ichan2")
        for seg in self.oriM1:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriM1.cm=1.4
        
        self.oriM2.insert("ichan2")
        for seg in self.oriM2:
            seg.gnatbar_ichan2 = gna  		# original 0.015
            seg.gkfbar_ichan2 = 0.013  		
            seg.gl_ichan2 = 0.00018
        self.oriM2.cm=1.4
        
        self.orit1.insert("ichan2")
        for seg in self.orit1:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
        self.orit1.cm=1.4
        
        self.orit2.insert("ichan2")
        for seg in self.orit2:
            seg.gnatbar_ichan2 = gna  		# Sodium conductance (original 0.015)
            seg.gkfbar_ichan2 = 0.013  		# Delayed K+ rectifier (fast)
            seg.gl_ichan2 = 0.00018         # Leak conductance
        self.orit2.cm=1.4
        
        for sec in self.all:
            # self.cm = Not setting cm
            
            sec.insert("ccanl")
            for seg in sec:
                seg.catau_ccanl = 10		# Time constant for decay of intracellular Ca2+
                seg.caiinf_ccanl = 5.e-6		# Steady-state intracellular Ca2+ concentration
            
            sec.insert("borgka")
            for seg in sec:
                seg.gkabar_borgka = 0.00015		# A-type K+ conductance

            sec.insert("nca") # N-type Ca2+ conductance
            for seg in sec:
                seg.gncabar_nca = 0.0008   		# check to modify- original 0.004
           
            sec.insert("lca") 
            for seg in sec:
                seg.glcabar_lca = 0.005		# L-type Ca2+ conductance
           
            sec.insert("gskch") 
            for seg in sec:
                seg.gskbar_gskch = 0.000002		# Ca2+-dependent K (SK) conductance
           
            sec.insert("mykca") 
            if (keepoldtypo==0):
                # commenting out gbar to match their typo
                for seg in sec:
                    seg.gkbar_mykca = 0.0002			# Ca2+ and Voltage-dependent K+ (BK) conductance
		 					# make catau slower70e-3 	cao=2 cai=50.e-6
            sec.Ra = 100			# 31.3 +/- 10.9
            sec.enat = 55
            sec.ekf = -90
            sec.ek = -90
            sec.elca = 130
            sec.esk = -90
            sec.el_ichan2 = -60			#-60.06         
        self.cao_ccanl = 2

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.nc.append(h.NetCon(self.soma(0.5)._ref_v, target, sec=self.soma))
        self.nc[-1].threshold = -10 # mV
        self.nc[-1].delay = delay # ms
        self.nc[-1].weight[0] = weight # NetCon weight is a vector    
        return self.nc[-1]

    def addSynapses(self):
        self.pre_list = []
       
        # E0
        syn_ = h.MyExp2Syn(self.radM1(0.5))
        self.pre_list.append(syn_)    # AMPA        EC (not used)
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E1
        syn_ = h.MyExp2Syn(self.radM2(0.5))
        self.pre_list.append(syn_)    # AMPA        EC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E2
        syn_ = h.MyExp2Syn(self.radM1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E3
        syn_ = h.MyExp2Syn(self.radM2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E4
        syn_ = h.MyExp2Syn(self.radT1(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E5
        syn_ = h.MyExp2Syn(self.radT2(0.5))
        self.pre_list.append(syn_)    # AMPA        CA3 Shaffer collateral
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E6
        syn_ = h.MyExp2Syn(self.oriT1(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0

        # E7
        syn_ = h.MyExp2Syn(self.oriT2(0.5))
        self.pre_list.append(syn_)    # AMPA        PC
        syn_.tau1 = 0.5
        syn_.tau2 = 3
        syn_.e = 0
        
        # I8
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A	Neighboring bistratified cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I9
        syn_ = h.MyExp2Syn(self.soma(0.5))
        self.pre_list.append(syn_)    # GABA-A	Basket cell
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I10
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I11
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-A	Septum
        syn_.tau1 = 1
        syn_.tau2 = 8
        syn_.e = -75
        
        # I12
        syn_ = h.MyExp2Syn(self.oriT1(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75
        
        # I13
        syn_ = h.MyExp2Syn(self.oriT2(0.6))
        self.pre_list.append(syn_)    # GABA-B	Septum
        syn_.tau1 = 35
        syn_.tau2 = 100
        syn_.e = -75

class BurstCell(modelcell):  
    """ Burst cell with stim attribute that references a BurstStim2 point process """
    def __init__(self, gid=-1):
        super().__init__()
        self.gid = gid

        self.is_art =1
        self.ncstim = []
        self.stim = h.BurstStim2()
       	self.stim.number = 10000
       	self.stim.start = 0
       	self.stim.interval = 10
       	self.stim.noise = 0
       	self.stim.burstint = 100	# interburst interval (ms)
       	self.stim.burstlen = 100	# burst length (ms)
           
    def __repr__(self):
        return "Burst cell {}: {} ms on, {} ms off, intraburst int is {} ms".format(self.gid, self.stim.burstint, self.stim.burstlen, self.stim.interval)

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.ncstim.append(h.NetCon(self.stim, target))
        self.ncstim[-1].delay = delay # ms
        self.ncstim[-1].weight[0] = weight # NetCon weight is a vector    
        return self.ncstim[-1]

class StimCell(modelcell):  
    """ Stim cell with stim attribute that references a RegnStim point process """
    def __init__(self, gid=-1):
        super().__init__()
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

    def connect2target(self, target, delay = 1, weight=0.04): # { localobj nc #$o1 target point process, optional $o2 returned NetCon
        self.ncstim.append(h.NetCon(self.stim, target))
        self.ncstim[-1].delay = delay # ms
        self.ncstim[-1].weight[0] = weight # NetCon weight is a vector    
        return self.ncstim[-1]