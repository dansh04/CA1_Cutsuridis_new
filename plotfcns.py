#################
# PLOT RESULTS
#################
import matplotlib.pyplot as plt
import fig9_patternrecall as fig9
import fig10_Vtraces as fig10
import numpy as np
import pickle

def plotresults(*args):
    """
    plotresults(simname)
    plotresults(params)
    
    Will generate plots of the simulation results

    Parameters
    ----------
    Optionally takes either of the following:
    simname: name of simulation
    params: parameter dictionary associated with the simulation

    Returns
    -------
    None.

    """
    
    for a in args:
        if type(a) is dict:
            params = a
        else:
            try:
                with open('pyresults/' + a + '.pickle', 'rb') as f:
                    params = pickle.load(f)
            except:
                print('*********** Warning! ***********')
                print('No parameter info found. Using default values:')
                params = {"simname":'par',"netfile":'N100S20P5',"numCycles":2,"network_scale":.2,"SIMDUR":550,"dt":0.025}
    for key in params:
        print(key,': ',params[key])
    
    fstem = "pyresults/" + params["simname"]
    spks = np.loadtxt("{}_spt.dat".format(fstem),skiprows=1)
    plt.figure()
    plt.scatter(spks[:,0],spks[:,1],s=.1)
    plt.xlabel("Time (ms)")
    plt.ylabel("Neuron #")
    plt.title("Spike Raster")
    plt.show()

    pvsoma = np.loadtxt("{}_pvsoma.dat".format(fstem),skiprows=1)
    plt.figure()
    plt.plot(pvsoma[:,0],pvsoma[:,1])
    plt.xlabel("Time (ms)")
    plt.ylabel("Membrane Potential (mV)")
    plt.title("Pattern Pyramidal Cell")
    plt.show()
        
    overall_performance=fig9.plot_results(params["simname"],params["netfile"],params["numCycles"],params["network_scale"])
    fig10.plot_voltages(params["simname"], 200, params["SIMDUR"] , params["dt"])
    print("overall_performance =",overall_performance)
    
    print( "** Finished plotting **")