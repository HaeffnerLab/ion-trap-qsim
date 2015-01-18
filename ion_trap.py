#!/usr/bin/env python
import numpy as np
from simulation_parameters import simulation_parameters
from equilibrium_positions import equilibrium_positions

class IonTrap:

        def __init__(self, omegax, omegaz, potential_config='harmonic'):
                
                self.omegax   =  2*np.pi * omegax
                self.omegaz   =  2*np.pi * omegaz
                self.potential_config = potential_config

        def load(self, chain):

                self.chain    =   chain
                #Find the eq positions:
                p             =   simulation_parameters()

                if self.potential_config == 'harmonic':
                        chain.set_zpositions( equilibrium_positions.get_positions(self.chain.num_of_ions, self.omegaz, p) )
                else:
                        print "Working on it."
                #Find the couplings:

                self.loaded = True
        
        def set_omegaz(self, omegaz):

                self.omegaz   =  2*np.pi * omegaz
                if self.loaded  == True:
                        self.load(self.chain)

        def set_omegax(self, omegax):
                
                self.omegax   =  2*np.pi * omegax
