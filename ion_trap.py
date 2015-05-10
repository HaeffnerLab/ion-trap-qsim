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
                if potential_config == 'harmonic':
                        self.potential   =  Potential(potential_config, self.omegaz)
                elif potential_config == 'generic':
                         self.potential   =  Potential(potential_config )

                #Frequency reference is the frequency with respect to
                # which all other frequencies in the entire simulation 
                # are measured:
                self.freq_reference  = freq_reference  
                self.magnetic_field_gradient = magnetic_field_gradient

        def load(self, chain, zpositions=[]):

                try:
                        potential_config = self.potential.config
                except AttributeError, e:
                        print("Potential is not set: IonTrap.set_potential()")

                #self.zpositions = zpositions
                self.chain    =   chain
                #Find the eq positions:
               

                #Set the position of ions:
                if self.potential.config == 'harmonic' or self.potential.config == 'generic':
                        self.chain.set_zpositions( equilibrium_positions.get_positions(self.chain.num_of_ions, self.potential) )
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

        def set_potential(self, fz ):
                """ Take fz = V(x_0, y_0, z), the potential along Z axis 
                    and pass it to the Potential instance.

                """

                self.potential.config = 'generic'
                self.potential.set_fz( fz )
               

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

                self.axial_freq  = sqrt(array( func )[::-1][2]/(0.5*self.p.mass))
                position_scale_factor  = (self.p.coulomb_coeff / (self.axial_freq**2 * self.p.mass))**(1./3) 
                print(position_scale_factor)
                self.Vz_deriv = polyder(func)
                print(self.Vz_deriv)
                print(self.Vz_deriv.order)
                self.Vz_deriv_inverted_arr = array(self.Vz_deriv)[::-1]
                print( self.Vz_deriv_inverted_arr  )
                self.rescaled_Vz_deriv_inverted_arr = array( [  self.Vz_deriv_inverted_arr[v]* (position_scale_factor**v) for v in range( len( self.Vz_deriv_inverted_arr ) ) ] )/(position_scale_factor*self.p.mass*self.axial_freq**2)
                
                #Constant term is important in rescaled_Vz_deriv since it represents electric field along Z axis and 
                #leads to a shift in trap center and positions along Z axis.
                #self.rescaled_Vz_deriv_inverted_arr[0] = 0. #Constant term is the last of the poly1d arrays
                self.rescaled_Vz_deriv_arr = self.rescaled_Vz_deriv_inverted_arr[::-1]
                self.rescaled_Vz_deriv = poly1d( self.rescaled_Vz_deriv_arr )
                print(self.axial_freq)

        def update(self, update_func):

                self.Vz += update_func

