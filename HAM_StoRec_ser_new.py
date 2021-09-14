# CA1 heteroassociative memory network: Storage and recall
# CA1 PCs, BCs, AACs, BSCs and OLMs (using moderate cell models)
# EC, CA3 (excitatory) and Septal (inhibitory) inputs
# Cycle is: Recall-Storage-Recall etc
# Serial code adapted from Hines' ran3ser.hoc
# VCU & BPG 8-1-09

# Results reported in V. Cutsuridis, S. Cobb and B.P. Graham,
# "Encoding and retrieval in a model of the hippocampal CA1 microcircuit",
# Hippocampus, in press, DOI 10.1002/hipo.20661, 2009.

#%%
#################
# MODULES
#################

# Import modules
from neuron import h, gui
import numpy as np
import random
import math
import sys
import netfcns
#%%

#################
# PARAMETERS
#################

# Set default values for parameters that can be passed in at the command line
plotflag = 1
batchflag = 1

scaleDown = 1
scaleEScon = 1

numCycles = 8
simname="test"

# Check for parameters being passed in via the command line
if len(sys.argv)>1:
    simname = sys.argv[1]
    if len(sys.argv)>2:
        batchflag = int(sys.argv[2])
        if len(sys.argv)>3:
            plotflag = int(sys.argv[3])
            if len(sys.argv)>4:
                scaleDown = float(sys.argv[4])
                if len(sys.argv)>5:
                    scaleEScon = float(sys.argv[5])
                    if len(sys.argv)>6:
                        numCycles = float(sys.argv[6])

fstem = "pyresults/" + simname
print("simname = ", simname, ", fstem = ", fstem)

# Set Timing Parameters
from model_const import *
from connprops import *

SPATT = calcSPATT(scaleDown)

h.tstop = SIMDUR
h.celsius = 34
#%%

#################
# CELL CREATION
#################

# Set GID ranges of cells and Load Cell Class definitions
class CellPop:
    def __init__(self, num=0, gidst=0, gidend=0, order=0, popname="", classtype="", isart=0):
        self.num=num
        self.gidst=gidst
        self.gidend=gidend
        self.order=order
        self.popname=popname
        self.classtype=classtype
        self.isart=isart

poplist=[]
poplist.append(CellPop(num=100*scaleDown, order=0, classtype="PyramidalCell", popname="PyramidalCell", isart=0))
poplist.append(CellPop(num=2, order=1, popname="BasketCell", classtype="BasketCell", isart=0))
poplist.append(CellPop(num=1, order=2, popname="AACell", classtype="AACell", isart=0))
poplist.append(CellPop(num=1, order=3, popname="BistratifiedCell", classtype="BistratifiedCell", isart=0))
poplist.append(CellPop(num=1, order=4, popname="OLMCell", classtype="OLMCell", isart=0))
poplist.append(CellPop(num=100*scaleDown, order=5, popname="CA3Cell", classtype="StimCell", isart=1))
poplist.append(CellPop(num=20*scaleDown, order=6, popname="ECCell", classtype="StimCell", isart=1))
poplist.append(CellPop(num=10*scaleDown, order=7, popname="SEPCell", classtype="BurstCell", isart=1))

dictpop={}

st=0
for pop in poplist:
    dictpop[pop.popname] = pop
    pop.gidst=st
    pop.gidend=pop.gidst + pop.num - 1
    st = pop.gidend + 1


from cellClasses import PyramidalCell, BasketCell,  \
AACell, BistratifiedCell, OLMCell, StimCell, BurstCell


# Calculate totals of cells
ncell = sum([o.num for o in poplist if o.isart==0])    # total number of real cells
nstim = sum([o.num for o in poplist if o.isart==1])        # total number of inputs
ntot = ncell+nstim # total number of cells


# creates the cells and appends them to a List called cells
# argument is the number of cells to be created
cells = []
ranlist = []
gidvec = []
nclist = []
h('{load_file("ranstream.hoc")}')  # to give each cell its own sequence generator

i=0
for pop in poplist:
    for j in range(int(pop.num)):    
        newcell = None
        exec("newcell = "+pop.classtype+ "(int("+str(i)+"))")
        cells.append(newcell)
        ranlist.append(h.RandomStream(i))  # ranlist.o(i) corresponds to
        gidvec.append(i)
        i+=1
#%%

#################
# SET STIMULATION
#################


# Configures the stimulation:

for i in range(int(dictpop["CA3Cell"].gidst),int(dictpop["CA3Cell"].gidend+1)):
    cells[i].stim.number = ENUM
    cells[i].stim.start = ESTART
    cells[i].stim.interval = EINT
    cells[i].stim.noise = ENOISE

for i in range(int(dictpop["ECCell"].gidst),int(dictpop["ECCell"].gidend+1)):
    cells[i].stim.number = ENUM
    cells[i].stim.start = ESTART
    cells[i].stim.interval = EINT
    cells[i].stim.noise = ENOISE


for i in range(int(dictpop["SEPCell"].gidst),int(dictpop["SEPCell"].gidend+1)):
    cells[i].stim.number = SEPNUM
    cells[i].stim.start = SEPSTART
    cells[i].stim.interval = SEPINT
    cells[i].stim.noise = SEPNOISE
    cells[i].stim.burstint = SEPBINT
    cells[i].stim.burstlen = SEPBLEN
    
    # Use the gid-specific random generator so random streams are
    # independent of where and how many stims there are.

#%%

#################
# SET STORAGE AND RECALL FUNCTIONALITY
#################

# Pattern storage and recall parameters

C_P = 1  # probability of excitatory connections received by each CA1 PC
         # from CA3 inputs (1 gives full connectivity)
         
SPATT = 20*scaleDown    # number of active cells per pattern
NPATT = 5    # number of stored patterns
NSTORE = 5    # number of new patterns to store

CPATT = 1    # index of cue pattern
CFRAC = 1    # fraction of active cells in cue

# file name of connection weights and patterns
# (cue and EC patterns taken from FSTORE file to implement storage)
# (use same file for FPATT and FSTORE to test recall only)
if scaleDown==1:
    FCONN = "weights/wgtsN100S20P5.dat" #"Weights/wgtsN100S20P5.dat"
    FPATT = "Weights/pattsN100S20P5.dat"    # "Weights/pattsN100S20P5.dat"    # already stored patterns
    FSTORE = "Weights/pattsN100S20P5.dat"    # "Weights/pattsN100S20P5.dat"    # new patterns to store
else:
    FCONN = "weights/wgtsN100S20P5scaled.dat" #"Weights/wgtsN100S20P5.dat"
    FPATT = "Weights/pattsN100S20P5scaled.dat"    # "Weights/pattsN100S20P5.dat"    # already stored patterns
    FSTORE = "Weights/pattsN100S20P5scaled.dat"    # "Weights/pattsN100S20P5.dat"    # new patterns to store
#%%

#################
# SET CONNECTIVITY
#################

# of connections received by each POSTCELLTYPE from a given PRECELLTYPE
class popConn:
    def __init__(self, popname="", prepop="", prenum=1, type="", weight=0, delay=3, synst=0, synend=0):
        self.popname=popname
        self.prepop=prepop
        self.prenum=prenum
        self.type=type
        self.weight=weight
        self.delay=delay
        self.synst=synst
        self.synend=synend
        
    def __repr__(self):
        return "Connection {} {} ---> {} of type {} with weight {}".format(int(self.prenum), self.prepop, self.popname, self.type, self.weight)

# Synapse indices
# onto CA1 PCs
E_EC = 0    # EC AMPA excit to medium SLM (2 of)
E_CA3 = 2    # CA3 AMPA excit to medium SR
EN_CA3 = 3    # CA3 NMDA excit to medium SR
EM_CA3 = 23    # CA3 modifiable (STDP) AMPA excit to medium SR
E_PC = 4    # CA1 recurrent AMPA excit to proximal SR
I_BC = 5    # ff&fb inhib via BCs to soma
I_AAC = 6    # ff&fb inhib via AACs to axon initial segment
I_BSC = 11    # ff&fb inhib via BSCs to SR med (12 of: 6 GABAA, 6 GABAB)
I_OLM = 7    # fb inhib via OLMs to SLM (4 of: 2 GABAA, 2 GABAB)

# onto INs (BC, AAC, BSC)
EI_EC = 0    # EC AMPA excit (2 of; not onto BSC)
EI_CA3 = 2    # CA3 AMPA excit (4 of)
EI_PC = 6    # CA1 PC AMPA excit (2 of)
II_SAME = 8    # inhib from neighbouring INs (BC->BC; BSC->BSC)
II_OPP = 9    # inhib from other INs (BSC->BC; BC->BSC)
II_SEP = 10    # inhib from septum (4 of: 2 GABAA, 2 GABAB)

# onto INs (OLM)
EO_PC = 0    # CA1 PC AMPA excit (2 of)
IO_IN = 2    # inhib from INs and septum (2 of: 1 GABAA, 1 GABAB)

connlist=[]                    # Postsynaptic type      Presynaptic type  number of connections (at least 1)
connlist.append(popConn(popname="PyramidalCell", prepop="CA3Cell", prenum=max([dictpop["CA3Cell"].num * scaleEScon, 1]), type="AMPA", weight=CLWGT, delay=CDEL, synst=E_CA3, synend=E_CA3+EN_CA3)) # CA3_PC. Use EM_CA3 for modifiable synapses
connlist.append(popConn(popname="BasketCell",    prepop="CA3Cell", prenum=max([dictpop["CA3Cell"].num * scaleEScon, 1]), type="AMPA", weight=EIWGT, delay=EIDEL, synst=EI_CA3, synend=EI_CA3+3)) # 
connlist.append(popConn(popname="AACell",        prepop="CA3Cell", prenum=max([dictpop["CA3Cell"].num * scaleEScon, 1]), type="AMPA", weight=EIWGT, delay=EIDEL, synst=EI_CA3, synend=EI_CA3+3)) # CA3_AAC
connlist.append(popConn(popname="BistratifiedCell", prepop="CA3Cell", prenum=max([dictpop["CA3Cell"].num * scaleEScon, 1]), type="AMPA", weight=EIWGT, delay=EIDEL, synst=EI_CA3, synend=EI_CA3+3)) # CA3_BSC

connlist.append(popConn(popname="PyramidalCell", prepop="ECCell", prenum=max([dictpop["ECCell"].num * scaleEScon, 1]), type="AMPA", weight=ECWGT, delay=EIDEL, synst=E_EC, synend=E_EC+2)) # EC_PC
connlist.append(popConn(popname="BasketCell",    prepop="ECCell", prenum=max([dictpop["ECCell"].num * scaleEScon, 1]), type="AMPA", weight=EIWGT, delay=EIDEL, synst=EI_EC, synend=EI_EC+1)) # EC_BC
connlist.append(popConn(popname="AACell",        prepop="ECCell", prenum=max([dictpop["ECCell"].num * scaleEScon, 1]), type="AMPA", weight=EIWGT, delay=EIDEL, synst=EI_EC, synend=EI_EC+1)) # EC_AAC

connlist.append(popConn(popname="BasketCell",    prepop="SEPCell", prenum=dictpop["SEPCell"].num, type="GABAA", weight=SEPWGT, delay=SEPDEL, synst=II_SEP, synend=II_SEP+1)) # SEP_BC
connlist.append(popConn(popname="AACell",        prepop="SEPCell", prenum=dictpop["SEPCell"].num, type="GABAA", weight=SEPWGT, delay=SEPDEL, synst=II_SEP, synend=II_SEP+1)) # SEP_AAC
connlist.append(popConn(popname="BistratifiedCell", prepop="SEPCell", prenum=dictpop["SEPCell"].num, type="GABAA", weight=SEPWGT, delay=SEPDEL, synst=II_SEP, synend=II_SEP+1)) # SEP_BSC
connlist.append(popConn(popname="OLMCell",       prepop="SEPCell", prenum=dictpop["SEPCell"].num, type="GABAA", weight=SEPWGT, delay=SEPDEL, synst=IO_IN, synend=IO_IN)) # SEP_OLM

connlist.append(popConn(popname="PyramidalCell", prepop="PyramidalCell", prenum=max([1 * scaleEScon, 1]), type="AMPA", weight=Pcell2Pcell_weight, delay=Pcell2Pcell_delay, synst=E_PC, synend=E_PC)) # PC_PC
connlist.append(popConn(popname="BasketCell",    prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = Pcell2Bcell_weight, delay = Pcell2Bcell_delay, synst=EI_PC, synend=EI_PC+1)) # PC_BC
connlist.append(popConn(popname="BistratifiedCell", prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = Pcell2AAcell_weight, delay = Pcell2AAcell_delay, synst=EI_PC, synend=EI_PC+1)) # PC_BSC
connlist.append(popConn(popname="AACell",        prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = Pcell2BScell_weight, delay = Pcell2BScell_delay, synst=EI_PC, synend=EI_PC+1)) # PC_AAC
connlist.append(popConn(popname="OLMCell",       prepop="PyramidalCell", prenum=max([dictpop['PyramidalCell'].num * scaleEScon, 1]), type="AMPA", weight = Pcell2OLMcell_weight, delay = Pcell2OLMcell_delay, synst=EO_PC, synend=EO_PC+1)) # PC_OLM

connlist.append(popConn(popname="PyramidalCell", prepop="BasketCell", prenum=2, type="GABAA", weight = Bcell2Pcell_weight, delay = Bcell2Pcell_delay, synst=I_BC, synend=I_BC)) # BC_PC
connlist.append(popConn(popname="BasketCell",    prepop="BasketCell", prenum=1, type="GABAA", weight = Bcell2Bcell_weight, delay = Bcell2Bcell_delay, synst=II_SAME, synend=II_SAME)) # BC_BC
connlist.append(popConn(popname="BistratifiedCell", prepop="BasketCell", prenum=2, type="GABAA", weight = Bcell2BScell_weight, delay = Bcell2BScell_delay, synst=II_OPP, synend=II_OPP)) # BC_BSC

# The call to make this connection is commented out in hoc code
# connlist.append(popConn(popname="OLMCell",       prepop="BasketCell", prenum=2, type="GABAA", weight = 0.0, delay = 1)) # BC_OLM

connlist.append(popConn(popname="PyramidalCell", prepop="AACell", prenum=1, type="GABAA", weight = AAcell2Pcell_weight, delay = AAcell2Pcell_delay, synst=I_AAC, synend=I_AAC)) # AAC_PC
connlist.append(popConn(popname="PyramidalCell", prepop="BistratifiedCell", prenum=1, type="GABAA", weight = BScell2Pcell_weight, delay = BScell2Pcell_delay, synst=I_BSC, synend=I_BSC+5)) # BSC_PC
connlist.append(popConn(popname="PyramidalCell", prepop="BistratifiedCell", prenum=1, type="GABAB", weight = BScell2Pcell_GABAB_weight, delay = BScell2Pcell_delay, synst=I_BSC+6, synend=I_BSC+11)) # BSC_PC

#connlist.append(popConn(popname="BistratifiedCell",    prepop="BistratifiedCell", prenum=1, type="GABAA", weight = BScell2BScell_weight, delay = BScell2BScell_delay, synst=II_SAME, synend=II_SAME)) # BSC_BC
connlist.append(popConn(popname="BasketCell",    prepop="BistratifiedCell", prenum=1, type="GABAA", weight = BScell2Bcell_weight, delay = BScell2Bcell_delay, synst=II_OPP, synend=II_OPP)) # BSC_BC

connlist.append(popConn(popname="PyramidalCell", prepop="OLMCell", prenum=1, type="GABAA", weight = OLMcell2Pcell_weight, delay = OLMcell2Pcell_delay, synst=I_OLM, synend=I_OLM+1)) # OLM_PC
connlist.append(popConn(popname="PyramidalCell", prepop="OLMCell", prenum=1, type="GABAB", weight =OLMcell2Pcell_GABAB_weight, delay = OLMcell2Pcell_delay, synst=I_OLM+2, synend=I_OLM+3)) # OLM_PC

# In the original Cutsuridis code, this weight is set first to 0.01 and then overwritten a few lines later to 0.0.
# The call to make this conection is commented out in hoc
# connlist.append(popConn(popname="BasketCell",    prepop="OLMCell", prenum=1, type="GABAA", weight = OLMcell2Bcell_weight, delay = OLMcell2Bcell_delay, synst=II_OPP, synend=II_OPP)) # OLM_BC
  
connect_random_low_start_ = 1  # low seed for mcell_ran4_init()
h.mcell_ran4_init(connect_random_low_start_)
nclist = []

# # Make connections with data from above
for conn in connlist: 
    print(conn)
    netfcns.connectcells(cells,ranlist, nclist, npost=dictpop[conn.popname].num, postgidstart=dictpop[conn.popname].gidst, npre = dictpop[conn.prepop].num, pregidstart = dictpop[conn.prepop].gidst, synstart=conn.synst, synend=conn.synend, npresyn=conn.prenum, weight=conn.weight, delay= conn.delay)




#netfcns.mkinputs(cells, dictpop['CA3Cell'].gidst, dictpop['ECCell'].gidst, dictpop['SEPCell'].gidst, ntot, dictpop)
# EC input to PCs
netfcns.connectEC(FPATT, ECPATT, NPATT, E_EC, 2, cells, dictpop['PyramidalCell'].gidst, dictpop['ECCell'].gidst, dictpop)	#  restore existing pattern
# CA3 input to PCs
netfcns.connectCA3(FCONN, C_P, EM_CA3, EN_CA3, cells, dictpop)	# with modifiable synapses
#%%

#################
# SET CUES FOR PATTERNS
#################

netfcns.mkcue(FPATT, CPATT, CFRAC, NPATT, SPATT, cells, dictpop['CA3Cell'].gidst, ranlist)	# cue from already stored pattern
#mkcue(FSTORE, CPATT, CFRAC, NSTORE)	# cue from new pattern
netfcns.mkEC(cells, ranlist, dictpop['ECCell'].gidst, dictpop["ECCell"].num)

#%%

#################
# DESIGN RESULTS
#################

netfcns.spikerecord(cells)
results = netfcns.vrecord(cells,dictpop, iPPC, iNPPC)



#%%

#################
# RUN SIMULATION AND WRITE RESULTS
#################

h('StepBy=100') # ms

h('walltime = startsw()')
h.xopen("midbalfcn.hoc")
h('objref fihw')
h('fihw = new FInitializeHandler(2, "midbal()")')



#if (batchflag==1):
print("Now running simulation at scale = ", scaleDown, " for time = ", h.SIMDUR, " with scaleEScon = ", scaleEScon)
h.run()
netfcns.spikeout(cells,fstem)
netfcns.vout(cells,results,fstem,dictpop)
print( "** Finished running sim and printing results **")

#%%

#################
# PLOT RESULTS
#################

if (plotflag==1):
    netfcns.spikeplot(cells,h.tstop,ntot)
    netfcns.vplot(cells,results)
 
print( "** Finished plotting **")
# else:
# 	# panel for simulation results
# 	h.xopen("HAM_SR.ses")
# 	xspikeres()




# At end
#spikeout(cells,fstem)
