# General framework - Test case 1: Giving an arbitrary position array:
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
from numpy import *
import numpy as np
from simulation_parameters import SimulationParameters, OptimizationParameters
from ions import Chain
from ion_trap import IonTrap




if __name__ == '__main__':

        N = 1
        chain     = Chain(N,2)
        trap      = IonTrap()
        trap.set_trap_radial_freq(1.8e6, 'Y')
        trap.set_trap_radial_freq(1.8e6, 'X')

        #Either set the electrode voltages
        trap.set_used_electrodes([4,5,6])
        trap.set_electrode_voltages([.1,-.2,.1])
        trap.load(chain) #--> Generate error if the potential cannot trap the given number of ions
                         # of if it leads to zigzag config
