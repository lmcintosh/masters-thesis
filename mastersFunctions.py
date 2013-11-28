# Master's Thesis functions

import numpy as np
from pylab import *
from scipy import *
from matplotlib import *
import pyentropy


# Adaptive exponential integrate-and-fire neuron
def aEIF(adaptationIndex, inputCurrent, v0):
    """
    Adaptive exponential integrate-and-fire neuron from Gerstner 2005
 
    Parameters
    ----------
    adaptationIndex : degree to which neuron adapts
                    
    inputCurrent : np array of the stimulus with size (M,N)
    
    v0 : specify the membrane potential at time 0
 
    Returns
    -------
    V : membrane voltage for the simulation
    w : adaptation variable for the simulation
    spikes : 0 or 1 for each time bin
    sptimes : list of spike times
    """
    
    # keep in mind that inputCurrent needs to start at 0 for sys to be in equilibrium
    
    # Physiologic neuron parameters from Gerstner et al.
    C       = 281    # capacitance in pF ... this is 281*10^(-12) F
    g_L     = 30     # leak conductance in nS
    E_L     = -70.6  # leak reversal potential in mV ... this is -0.0706 V
    delta_T = 2      # slope factor in mV
    V_T     = -50.4  # spike threshold in mV
    tau_w   = 144    # adaptation time constant in ms
    V_peak  = 20     # when to call action potential in mV
    b       = 0.0805 # spike-triggered adaptation
    a       = adaptationIndex
    noiseT  = 0.3
    
    # Simulation parameters
    delta = 0.5                      # dt
    M     = inputCurrent.shape[0]    # number of neurons
    N     = inputCurrent.shape[1]    # number of simulation points is determined by size of inputCurrent
    T     = np.linspace(0,N*delta,N) # time points corresponding to inputCurrent (same size as V, w, I)
    
    # Thermodynamic parameters
    kB   = 1.3806503*10**(-23)   # Boltzmann's constant
    beta = 1/(kB*310.65)        # kB times T where T is in Kelvin
    
    # Initialize variables
    V       = np.zeros((N,M))
    w       = np.zeros((N,M))
    spikes  = np.zeros((N,M))
    sptimes = [[] for _ in xrange(M)]
    V[0]    = v0                # this gives us a chance to say what the membrane voltage starts at
                                # so we can draw initial conditions from the Boltzmann dist. later
    
    # Run model
    for i in xrange(N-1):
        V[i+1,:] = V[i,:] + (delta/C)*( -g_L*(V[i,:] - E_L) + g_L*delta_T*np.exp((V[i,:] - V_T)/delta_T) - w[i,:] + inputCurrent[:,i+1]) + np.sqrt(noiseT)*randn(1,M)
        w[i+1,:] = w[i,:] + (delta/tau_w)*(a*(V[i,:] - E_L) - w[i,:])
    
        # spiking mechanism
        ind = where(V[i+1,:] >= V_peak)
        if size(ind[0]) > 0:
            V[i+1,ind]      = E_L
            w[i+1,ind]      = w[i,ind] + b
            spikes[i+1,ind] = 1
            [sptimes[j].append(T[i+1]) for j in ind[0]]
    
    return [V.transpose(),w.transpose(),spikes.transpose(),sptimes]    


# Define raster plot
def raster(event_times_list, color='k'):
    """
    Creates a raster plot
 
    Parameters
    ----------
    event_times_list : iterable
                       a list of event time iterables
    color : string
            color of vlines
 
    Returns
    -------
    ax : an axis containing the raster plot
    """
    ax = gca()
    for ith, trial in enumerate(event_times_list):
        vlines(trial, ith + .5, ith + 1.5, color=color,linewidth=1.3)
    ylim(.5, len(event_times_list) + .5)
    
    return ax



# Define input current types
def stepFn(meanCurrent,variance,N, M):
    ''' stepFn(meanCurrent, variance, N,M) '''
    
    x    = meanCurrent*np.ones((N,M)) + np.sqrt(variance)*np.random.randn(N,M)
    x[0] = 0
    
    return x.transpose()


def whiteNoise(meanCurrent,variance,N,M):
    ''' whiteNoise(meanCurrent, variance, N,M) '''
    
    x    = meanCurrent + np.sqrt(variance)*np.random.randn(N,M)
    x[0] = 0
    
    return x.transpose()


def shotNoise(meanCurrent,variance,N,M):
    ''' shotNoise(meanCurrent, variance, N,M) '''
    
    dt   = 0.5
    tau  = 3 # ms
    F    = 30
    t    = linspace(0,F*dt,F)
    x    = meanCurrent + np.sqrt(variance)*np.random.randn(N+F-1,M)
    y    = t*np.exp(-t/tau)
    z    = np.zeros((N,M))
    for i in xrange(M):
        z_i    = np.convolve(x[:,i], y, mode = 'valid')
        z[:,i] = z_i
    z    = z - np.mean(z) + meanCurrent
    z[0] = 0
    
    return z.transpose()


def ouProcess(meanCurrent,variance,N,M):
    ''' shotNoise(meanCurrent, variance, N,M) '''
    
    dt    = 0.5
    tau   = 3 # ms
    tau2  = N/20
    F     = 30
    t     = linspace(0,F*dt,F)
    t2    = linspace(0,N*dt,N)
    x     = meanCurrent + np.sqrt(variance)*np.random.randn(N+F-1,M)
    y     = np.exp(-t/tau)
    y2    = np.exp(-t2/tau2)
    z0    = 0
    drift = (z0 - meanCurrent)*y2
    z     = np.zeros((N,M))
    for i in xrange(M):
        z_i    = np.convolve(x[:,i], y, mode = 'valid')
        z[:,i] = z_i - np.mean(z_i) + meanCurrent + drift
    
    return z.transpose()

"""
brownian() implements one dimensional Brownian motion (i.e. the Wiener process).
"""
#from math import sqrt
from scipy.stats import norm
#import numpy as np


def brownian(x0, n, dt, delta, out=None):
    """\
    Generate an instance of Brownian motion (i.e. the Wiener process):

        X(t) = X(0) + N(0, delta**2 * t; 0, t)

    where N(a,b; t0, t1) is a normally distributed random variable with mean a and
    variance b.  The parameters t0 and t1 make explicit the statistical
    independence of N on different time intervals; that is, if [t0, t1) and
    [t2, t3) are disjoint intervals, then N(a, b; t0, t1) and N(a, b; t2, t3)
    are independent.
    
    Written as an iteration scheme,

        X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)


    If `x0` is an array (or array-like), each value in `x0` is treated as
    an initial condition, and the value returned is a numpy array with one
    more dimension than `x0`.

    Arguments
    ---------
    x0 : float or numpy array (or something that can be converted to a numpy array
         using numpy.asarray(x0)).
        The initial condition(s) (i.e. position(s)) of the Brownian motion.
    n : int
        The number of steps to take.
    dt : float
        The time step.
    delta : float
        delta determines the "speed" of the Brownian motion.  The random variable
        of the position at time t, X(t), has a normal distribution whose mean is
        the position at time t=0 and whose variance is delta**2*t.
    out : numpy array or None
        If `out` is not None, it specifies the array in which to put the
        result.  If `out` is None, a new numpy array is created and returned.

    Returns
    -------
    A numpy array of floats with shape `x0.shape + (n,)`.
    
    Note that the initial value `x0` is not included in the returned array.
    """

    x0 = np.asarray(x0)

    # For each element of x0, generate a sample of n numbers from a
    # normal distribution.
    r = norm.rvs(size=x0.shape + (n,), scale=delta*sqrt(dt))

    # If `out` was not given, create an output array.
    if out is None:
        out = np.empty(r.shape)

    # This computes the Brownian motion by forming the cumulative sum of
    # the random samples. 
    np.cumsum(r, axis=-1, out=out)

    # Add the initial condition.
    out += np.expand_dims(x0, axis=-1)

    return out



def ensemble(adaptiveIndex, numNeurons, inputType, duration):
    ''' Simulates an ensemble of neurons
    
    inputs: adaptiveIndex, numNeurons, and inputType
    
    possible inputType includes:
        'step'
        'white'
        'ou'
        'shotnoise'
        'brownian'
    
    returns: V (neurons x time), w, spikes, stimulus
    '''
    
    # constants
    N           = duration
    a           = adaptiveIndex
    M           = numNeurons
    v0          = -70.6    # mV
    meanCurrent = 555.0    # pA ?
    variance    = 100.0
    delta       = 0.5
    x0          = 0.0
    
    # initialize variables
    T = np.linspace(0,N*delta,N) # time points corresponding to inputCurrent (same size as V, w, I)
    
    if inputType == 'step':
        current = stepFn(meanCurrent,variance,N,M)
    elif inputType == 'white':
        current = whiteNoise(meanCurrent,variance,N,M)
    elif inputType == 'ou':
        current = ouProcess(meanCurrent,variance,N,M)
    elif inputType == 'shotnoise':
        current = shotNoise(meanCurrent,variance,N,M)
    elif inputType == 'brownian':
        current = brownian(np.ones((M,)) * meanCurrent, N, delta, variance/25.0, out=None)
        current[:,0] = x0
        
    V,w,spikes,sptimes = aEIF(a, current, v0)
        
    return V, w, spikes, sptimes, T, current



def mutiN(voltage, current, nBins, minV, maxV, minC, maxC):
    '''Function to compute the mutual information between neuron voltage and the input current - with just two bins!
    
    Note that you need to use asarray(V) and asarray(current) to get time slices instead of per-neuron slices'''
    
    # could potentially ask for min, max, and then use linspace to specify bin edges
    binsA = linspace(minV, maxV, nBins+1)
    binsB = linspace(minC, maxC, nBins+1)
    
    # note that in this function's MATLAB counterpart, I explicitly made one bin spiking.
    # here I avoid doing that because I don't force a spike to have a particular voltage,
    # I just call it when it passes threshold on the next iteration
    
    # Counting
    CountsAB, xedges, yedges = histogram2d(squeeze(voltage),squeeze(current), bins=[binsA, binsB])
    CountsA = sum(CountsAB,0) # this sums all the rows, leaving the marginal distribution of voltage
    CountsB = sum(CountsAB,1) # this sums all the cols, leaving the marginal distribution of current
    
    pA  = CountsA.astype(float)/float(sum(CountsA))
    pB  = CountsB.astype(float)/float(sum(CountsB))
    pAB = CountsAB.astype(float)/float(sum(CountsAB))
    
    
    if abs(1 - sum(pA)) > 0.01:
        print 'Probabilities pA do not sum to one ' + str(sum(pA))
    if abs(1 - sum(pB)) > 0.01:
        print 'Probabilities pB do not sum to one ' + str(sum(pB))
    if abs(1 - sum(pAB)) > 0.01:
        print 'Probabilities pAB do not sum to one ' + str(sum(pAB))
    
    
    # Entropies
    HA  = 0.0
    HB  = 0.0
    HAB = 0.0
    
    for prob in pB:
        if prob != 0:
            HB = HB - prob*log2(prob)
    
    for j in xrange(len(pA)):
        if pA[j] != 0:
            HA = HA - pA[j] * log2(pA[j])
            for i in xrange(len(pB)):
                if pAB[i,j] != 0:
                    HAB = HAB - pAB[i,j] * log2(pAB[i,j])
                    
    H = [HA, HB, HAB]
    I = HA + HB - HAB
    
    return H, I



def binaryWordsInformation(spikes,stimulus):
    '''Compute entropy of spike trains with binary words approach.
    
    Spikes and stimulus are both 1-d vertical numpy arrays
    with as many elements as neurons.
    '''
    # Quantize stimulus and
    # align spikes & stim for entropy calculation
    numNeurons = len(spikes)
    spikes     = spikes.flatten()   # get in format (numNeurons,)
    stimulus   = stimulus.flatten() # get in format (numNeurons,)
    if len(np.unique(spikes)) > 2:
        spikes = pyentropy.quantise(spikes, 2, uniform='sampling', minmax=None, centers=True)
        spikes = spikes[0]
    else:
        spikes     = spikes.astype(int) # make integer array from float64 array
    stim_quant = pyentropy.quantise(stimulus, 2, uniform='sampling', minmax=None, centers=True)
    system     = pyentropy.DiscreteSystem(spikes, (1,numNeurons), stim_quant[0], (1,numNeurons))
    
    # compute entropies
    system.calculate_entropies(method='plugin', calc=['HX', 'HXY'])
    
    # return information between spikes and stimulus
    return system.I()
