#A collection of methods to find :

import numpy as np 
import matplotlib.pylab as plt
import simulation_parameters
from numpy import absolute, real, linspace
from qutip import *

class Dynamics:
    
    def __init__(self, chain, time_precision):
        self.chain = chain
        self.time_precision = time_precision


    def construct_propagator(self, H, times, flag='all'):
        """Retrun propagator for Hamiltonian H. If flag=='all' return
        propagator at all time points given in times, if flag=='final' return
        propagator at time zero and times[-1].

        """

        if flag == 'all':
            time_evol_op_list       =  []
            time_evol_op            =  (-1.j*H* self.time_precision/10.).expm() #tensor(qeye(2), qeye(2), qeye(M), qeye(M))
            exponential             =  (-1.j*H* self.time_precision ).expm() 
            for t in enumerate(times):
                time_evol_op_list      +=   [time_evol_op]
                time_evol_op            =   exponential * time_evol_op

        elif flag == 'final':
            time_evol_op            =  identity_op #nntensor(qeye(2), qeye(2), qeye(M), qeye(M))
            exponential             =  (-1.j*H* self.time_precision ).expm() 
            for t in enumerate(times):
                time_evol_op            =   exponential * time_evol_op
            time_evol_op_list       =   [time_evol_op]
        
        self.time_evol_op_list  =  time_evol_op_list
            

