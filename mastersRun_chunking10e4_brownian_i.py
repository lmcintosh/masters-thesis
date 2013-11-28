# Master's Thesis iPython Notebook scripts to run

import numpy as np
from pylab import *
from scipy import *
from matplotlib import *
import datetime
import pickle

from mastersFunctions import *
from emailing import *


try:
    stim_type = 'brownian'
    version   = 0

    M     = int(10e4)
    nBins = 250
    adapt = np.linspace(0,25,100)
    if version == 0:
        adapt = adapt[0:30]
    elif version == 1:
        adapt = adapt[30:60]
    elif version == 2:
        adapt = adapt[60:100]
    
    steady_state_mem   = []
    steady_state_mem2  = []
    steady_state_pred  = []
    steady_state_pred2 = []
    steady_state_max   = []
    steady_state_max2  = []
    sum_mem            = []
    sum_mem2           = []
    sum_pred           = []
    sum_pred2          = []
    sum_max            = []
    sum_max2           = []
    i_mem_fnA          = np.zeros((3,999))
    i_mem_fnA2         = np.zeros((3,999))
    i_pred_fnA         = np.zeros((3,999))
    i_pred_fnA2        = np.zeros((3,999))
    i_max_fnA          = np.zeros((3,999))
    i_max_fnA2         = np.zeros((3,999))
    
    for a in adapt:
        # run simulation
        V,w,spikes,sptimes,T,stimulus = ensemble(a,int(M),stim_type,1000)
        
        V        = np.asarray(V)
        stimulus = np.asarray(stimulus)
        
        # compute i_mem, i_pred, i_max
        i_mem   = []
        i_mem2  = []
        i_pred  = []
        i_pred2 = []
        i_max   = []
        i_max2  = []

        # I used to flatten V, stimulus but this not only takes a long time
        # but gives unfair advantage to nonspiking Vs when the range is smaller
        # and discrimination is effectively much higher
        minV = -73.0
        maxV = 20.0
        minC = 0.0
        maxC = 880.0

        for i in xrange(shape(V)[1]-1):
            H, I = mutiN(V[:,i],stimulus[:,i],nBins,minV,maxV,minC,maxC)
            I2   = binaryWordsInformation(spikes[:,i],stimulus[:,i])
            i_mem.append(I)
            i_mem2.append(I2)
            H, I = mutiN(V[:,i],stimulus[:,i+1],nBins,minV,maxV,minC,maxC)
            I2   = binaryWordsInformation(spikes[:,i],stimulus[:,i+1])
            i_pred.append(I)
            i_pred2.append(I2)
            H, I = mutiN(stimulus[:,i],stimulus[:,i+1],nBins,minC,maxC,minC,maxC)
            I2   = binaryWordsInformation(stimulus[:,i],stimulus[:,i+1])
            i_max.append(I)
            i_max2.append(I2)
    
        steady_state_mem.append(i_mem[-1])
        steady_state_pred.append(i_pred[-1])
        steady_state_max.append(i_max[-1])
        steady_state_mem2.append(i_mem2[-1])
        steady_state_pred2.append(i_pred2[-1])
        steady_state_max2.append(i_max2[-1])
        
        sum_mem.append(sum(i_mem))
        sum_pred.append(sum(i_pred))
        sum_max.append(sum(i_max))

        sum_mem2.append(sum(i_mem2))
        sum_pred2.append(sum(i_pred2))
        sum_max2.append(sum(i_max2))
        
        if a < 0.1:
            i_mem_fnA[0,:]  = asarray(i_mem)
            i_pred_fnA[0,:] = asarray(i_pred)
            i_max_fnA[0,:]  = asarray(i_max)
            
            i_mem_fnA2[0,:]  = asarray(i_mem2)
            i_pred_fnA2[0,:] = asarray(i_pred2)
            i_max_fnA2[0,:]  = asarray(i_max2)
        elif abs(a - 6) < 0.1:
            i_mem_fnA[1,:]  = asarray(i_mem)
            i_pred_fnA[1,:] = asarray(i_pred)
            i_max_fnA[1,:]  = asarray(i_max)

            i_mem_fnA2[1,:]  = asarray(i_mem2)
            i_pred_fnA2[1,:] = asarray(i_pred2)
            i_max_fnA2[1,:]  = asarray(i_max2)
        elif abs(a - 25) < 0.1:
            i_mem_fnA[2,:]  = asarray(i_mem)
            i_pred_fnA[2,:] = asarray(i_pred)
            i_max_fnA[2,:]  = asarray(i_max)

            i_mem_fnA2[2,:]  = asarray(i_mem2)
            i_pred_fnA2[2,:] = asarray(i_pred2)
            i_max_fnA2[2,:]  = asarray(i_max2)
    
    

        with open('masters_data_steadyState_' + str(stim_type) + '_a' + str(a) + '_' + str(datetime.date.today()) + '.pik', 'wb') as f:
            pickle.dump([steady_state_mem, steady_state_pred, steady_state_max], f, -1)

        with open('masters_data_sum_' + str(stim_type) + '_a' + str(a) + '_' + str(datetime.date.today()) + '.pik', 'wb') as f:
            pickle.dump([sum_mem, sum_pred, sum_max], f, -1)

        with open('masters_data_examples_' + str(stim_type) + '_a' + str(a) + '_' + str(datetime.date.today()) + '.pik', 'wb') as f:
            pickle.dump([i_mem_fnA, i_pred_fnA, i_max_fnA], f, -1)

        with open('masters_data_steadyState2_' + str(stim_type) + '_a' + str(a) + '_' + str(datetime.date.today()) + '.pik', 'wb') as f:
            pickle.dump([steady_state_mem2, steady_state_pred2, steady_state_max2], f, -1)

        with open('masters_data_sum2_' + str(stim_type) + '_a' + str(a) + '_' + str(datetime.date.today()) + '.pik', 'wb') as f:
            pickle.dump([sum_mem2, sum_pred2, sum_max2], f, -1)

        with open('masters_data_examples2_' + str(stim_type) + '_a' + str(a) + '_' + str(datetime.date.today()) + '.pik', 'wb') as f:
            pickle.dump([i_mem_fnA2, i_pred_fnA2, i_max_fnA2], f, -1)



    emailWhenDone()
    
except:
    import traceback
    import smtplib
    
    sysexecinfo = sys.exc_info()
    
    emailWhenError(sysexecinfo)




