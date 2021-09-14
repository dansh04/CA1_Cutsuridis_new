# -*- coding: utf-8 -*-
"""
Created on Sat Aug  1 13:05:11 2020

@author: mbezaire
"""

STARTDEL = 50    # msecs
THETA = 250    # msecs (4 Hz)
GAMMA = 25    # msecs (40 Hz)
ECCA3DEL = 9    # msecs

# Septal inhibition
SEPNUM = 1000    # number of SEP spikes
SEPSTART = STARTDEL+(THETA/12)    # time of first SEP spike
SEPINT = 20    # SEP spike ISI (during burst)
SEPNOISE = 0.4    # SEP ISI noise
SEPBINT = 2*THETA/3    # SEP interburst interval
SEPBLEN = THETA/3    # SEP burst length
SEPWGT = 0.02	# SEP weight to BCs and AACs
SEPWGTL = 0.0002	# SEP weight to BSCs and OLMs
SEPDEL = 1	# SEP delay

# Background excitation (not used)
ENUM = 0    # number of spikes
ESTART = 0    # time of first spike
EINT = 200    # spike ISI
ENOISE = 1    # ISI noise
EWGT = 0.001    # excitatory weights (AMPA)
ENWGT = 0.002    # excitatory weights (NMDA)
EDEL = 1    # delay (msecs)

# EC excitation
ECPATT = 1    # index of output pattern
ECNUM = 1000    # number of EC spikes
ECSTART = STARTDEL    # time of first EC spike
ECINT = GAMMA    # EC spike ISI
ECNOISE = 0.2    # EC ISI noise
ECWGT = 0.0    # EC weight to PCs
#ECWGT = 0.001    # EC weight to PCs
ECDEL = 1    # EC delay
EIWGT = 0.00015    # excitatory weights to INs
EIDEL = 1    # delay (msecs)

# Cue (CA3) excitation
CNUM = 1000    # number of cue spikes
CSTART = STARTDEL+ECCA3DEL    # time of first cue spike
CINT = GAMMA    # cue spike ISI
CNOISE = 0.2    # cue ISI noise
CHWGT = 0.0015    # cue weight
CLWGT = 0.0005    # unlearnt weight (usually 0)
CNWGT = 0.0005    # excitatory weights (NMDA)
CDEL = 1    # cue delay

# STDP configuration
STDPDFAC = 0    # depression factor
STDPPFAC = 0    # potentiation factor
#STDPDFAC = 0.2    # depression factor
#STDPPFAC = 1.0    # potentiation factor
AMPASUPP = 0.4    # fraction of AMPA during storage
STDPTHRESH = -55    # voltage threshold for STDP
STDPSTART = STARTDEL+(THETA/2)    # STDP starts at same time as EC input
STDPINT = THETA/2    # STDP interburst (recall) interval
STDPLEN = THETA/2    # STDP burst (storage) length

C_P = 1  # probability of excitatory connections received by each CA1 PC
         # from CA3 inputs (1 gives full connectivity)
def calcSPATT(scaleDown):        
    SPATT = 20*scaleDown	# number of active cells per pattern
    return SPATT

NPATT = 5	# number of stored patterns
NSTORE = 5	# number of new patterns to store

CPATT = 1	# index of cue pattern
CFRAC = 1	# fraction of active cells in cue
iPPC=1		# index of a pattern PC (1st patt in 5 patterns)
iNPPC=0		# index of a non-pattern PC (1st patt in 5 patterns)

# Septal inhibition
SEPNUM = 1000    # number of SEP spikes
SEPSTART = STARTDEL+(THETA/12)    # time of first SEP spike
SEPINT = 20    # SEP spike ISI (during burst)
SEPNOISE = 0.4    # SEP ISI noise
SEPBINT = 2*THETA/3    # SEP interburst interval
SEPBLEN = THETA/3    # SEP burst length
SEPWGT = 0.02	# SEP weight to BCs and AACs
SEPWGTL = 0.0002	# SEP weight to BSCs and OLMs
SEPDEL = 1	# SEP delay

# Background excitation (not used)
ENUM = 0    # number of spikes
ESTART = 0    # time of first spike
EINT = 200    # spike ISI
ENOISE = 1    # ISI noise
EWGT = 0.001    # excitatory weights (AMPA)
ENWGT = 0.002    # excitatory weights (NMDA)
EDEL = 1    # delay (msecs)

# EC excitation
ECPATT = 1    # index of output pattern
ECNUM = 1000    # number of EC spikes
ECSTART = STARTDEL    # time of first EC spike
ECINT = GAMMA    # EC spike ISI
ECNOISE = 0.2    # EC ISI noise
ECWGT = 0.0    # EC weight to PCs
#ECWGT = 0.001    # EC weight to PCs
ECDEL = 1    # EC delay
EIWGT = 0.00015    # excitatory weights to INs
EIDEL = 1    # delay (msecs)

# Cue (CA3) excitation
CNUM = 1000    # number of cue spikes
CSTART = STARTDEL+ECCA3DEL    # time of first cue spike
CINT = GAMMA    # cue spike ISI
CNOISE = 0.2    # cue ISI noise
CHWGT = 0.0015    # cue weight
CLWGT = 0.0005    # unlearnt weight (usually 0)
CNWGT = 0.0005    # excitatory weights (NMDA)
CDEL = 1    # cue delay

# STDP configuration
STDPDFAC = 0    # depression factor
STDPPFAC = 0    # potentiation factor
#STDPDFAC = 0.2    # depression factor
#STDPPFAC = 1.0    # potentiation factor
AMPASUPP = 0.4    # fraction of AMPA during storage
STDPTHRESH = -55    # voltage threshold for STDP
STDPSTART = STARTDEL+(THETA/2)    # STDP starts at same time as EC input
STDPINT = THETA/2    # STDP interburst (recall) interval
STDPLEN = THETA/2    # STDP burst (storage) length
