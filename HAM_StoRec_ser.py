# CA1 heteroassociative memory network: Storage and recall
# CA1 PCs, BCs, AACs, BSCs and OLMs (using moderate cell models)
# EC, CA3 (excitatory) and Septal (inhibitory) inputs
# Cycle is: Recall-Storage-Recall etc
# Serial code adapted from Hines' ran3ser.hoc
# VCU & BPG 8-1-09

# Results reported in V. Cutsuridis, S. Cobb and B.P. Graham,
# "Encoding and retrieval in a model of the hippocampal CA1 microcircuit",
# Hippocampus, in press, DOI 10.1002/hipo.20661, 2009.

from neuron import h, gui
import numpy as np
import random
import math
import sys

h("strdef simname")
h("batchflag = 0")
h("plotflag = 0")
h("scaleDown = 1")
h("scaleEScon = 1") # Scaling argument for calcium channel conductances

h.simname="test"
if len(sys.argv)>1:
    h.simname = sys.argv[1]
    if len(sys.argv)>2:
        h.batchflag = int(sys.argv[2])
        if len(sys.argv)>3:
            h.plotflag = int(sys.argv[3])
            if len(sys.argv)>4:
                h.scaleDown = float(sys.argv[4])
                if len(sys.argv)>5:
                    h.scaleEScon = float(sys.argv[5])

print("h.simname = ", h.simname)


h.xopen("HAM_StoRec_ser_new.hoc")