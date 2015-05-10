#!/usr/bin/env python
from __future__ import division, absolute_import, print_function, unicode_literals
import numpy as np
from numpy import *
from simulation_parameters import simulation_parameters
from equilibrium_positions import equilibrium_positions


class IonTrap:

        def __init__(self, omegax, omegaz, potential_config='harmonic', freq_reference = 0.0, magnetic_field_gradient = 0.0):
                
                self.omegax   =  2*np.pi * omegax
                self.omegaz   =  2*np.pi * omegaz
                self.p             =   simulation_parameters()
                self.potential   =  Potential(potential_config, self.omegaz)

                #Frequency reference is the frequency with respect to
                # which all other frequencies in the entire simulation 
                # are measured:
                self.freq_reference  = freq_reference  
                self.magnetic_field_gradient = magnetic_field_gradient

        def load(self, chain, zpositions=[]):

                #self.zpositions = zpositions
                self.chain    =   chain
                #Find the eq positions:
               

                #Set the position of ions:
                if self.potential.config == 'harmonic' or self.potential.config == 'generic':
                        self.chain.set_zpositions( equilibrium_positions.get_positions(self.chain.num_of_ions, self.potential, self.p) )
                elif self.potential.config == 'positions':  
                        if len(zpositions) == self.chain.num_of_ions:
                                self.chain.set_zpositions( zpositions )
                        else:
                                raise Exception("Z position of ions must be given as a list with length equal to\n number of ions")

                else:
                        raise Exception('Wrong specification for potential_config {}.'.format(self.potential_config))
               
                self.chain.set_carrier_frequencies(self.freq_reference, self.magnetic_field_gradient)

                #Find and set the ion-ion couplings:
                self.chain.set_couplings( self.omegax )

                #Set expansion of chain's motion local in terms of normal destruction operators and vice versa, and eigenvalues.
                self.chain.set_normal_mode_structure()


                self.loaded   =  True
        
        def set_omegaz(self, omegaz):

                self.omegaz   =  2*np.pi * omegaz
                if self.loaded:
                        self.load(self.chain)

        def set_omegax(self, omegax):
                
                self.omegax   =  2*np.pi * omegax

        def set_ions_positions(self, zpositions):
                self.potential.config = 'positions'
                self.load(self.chain, zpositions)

        def set_potential(self, ):
                self.potential.config = 'generic'


class Potential(object):

        def __init__(self, potential_config, axial_freq=0):
                """ Contains information about total potential which is devided by 
                    0.5*M*omegaz**2 along the Z axis. 
                    func is a numpy.polynomial object, determining V(x_0,y_0, z), potential along the Z trap axis

                """ 
                self.p             =   simulation_parameters()

                self.config = potential_config

                if self.config == 'harmonic':
                        if axial_freq == 0:
                                raise Exception("Axial frequency is not set.")
                        self.axial_freq = axial_freq #Assuming that it is already multiplied by 2*np.pi before passing to Potential.
                
                
        def set_fz(self, func):
                """ Set func which is a numpy.polynomial (poly1d) object, determining V(x_0,y_0, z), potential along the Z trap axis
                """
                self.config = 'generic'
                self.axial_freq  = sqrt(array( self.func )[::-1][2]/(0.5*self.p.mass))
                position_scale_factor  = (self.p.coulomb_coeff / (omegaz**2 * self.p.mass))**(1./3) 
                self.Vz_deriv = polyder(self.func)
                self.Vz_deriv_arr = array(self.Vz_deriv)
                self.Vz_deriv_rescaled_arr = array( [  self.Vz_deriv_arr[v]* (position_scale_factor**v) for v in range( len( self.Vz_deriv_arr ) ) ] )/(position_scale_factor*self.p.mass*self.axial_freq**2)

                #Getting rid of constant term in case it is large (not good for minimization):
                self.Vz_deriv_rescaled_arr[-1] = 0 #Constant term is the last of the poly1d arrays
                self.Vz_deriv_rescaled = poly1d( self.Vz_deriv_rescaled_arr )

        def update(self, update_func):

                self.Vz += update_func
