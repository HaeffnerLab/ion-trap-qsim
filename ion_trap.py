#!/usr/bin/env python
#from __future__ import division, absolute_import, print_function, unicode_literals
import numpy as np
from numpy import *
from simulation_parameters import SimulationParameters
from equilibrium_positions import equilibrium_positions
from simulation_parameters import OptimizationParameters


class IonTrap:

        def __init__(self, freq_reference = 0.0, magnetic_field_gradient = 0.0, laser_axis='Y'):
                
                self.p             =   SimulationParameters()

                self.potential     =   Potential('generic')
                
                self.laser_axis    =  laser_axis
                #Frequency reference is the frequency with respect to
                # which all other frequencies in the entire simulation 
                # are measured:
                self.freq_reference  = freq_reference  
                self.magnetic_field_gradient = magnetic_field_gradient


        def load(self, chain):

                self.chain    =   chain
                

                #Set the position of ions:
                if self.potential.config == 'harmonic' or self.potential.config == 'generic':
                        self.chain.set_zpositions( equilibrium_positions.get_positions(self.chain.num_of_ions, self.potential) )
                
                self.chain.set_carrier_frequencies(self.freq_reference, self.magnetic_field_gradient)

                #Find and set the ion-ion couplings:
                self.chain.set_couplings( self.potential.get_trap_radial_freq(self.laser_axis) )

                #Set expansion of chain's motion local in terms of normal destruction operators and vice versa, and eigenvalues.
                self.chain.set_normal_mode_structure()


                self.loaded   =  True
        
        
        def reload(self):

                self.load(self.chain)


        def set_ions_positions(self, chain, zpositions):
                
                self.potential.config = 'positions'
                self.chain = chain

                if len(zpositions) == self.chain.num_of_ions:
                        self.chain.set_zpositions( zpositions )
                else:
                        raise Exception("Z position of ions must be given as a list with length equal to\n number of ions")


        '''
        def set_ions_positions(self):

                self.params.dummy_trap.set_potential(self.potential) #Set trap potential
                self.params.dummy_trap.load( self.params.chain ) #load ions with this potential
                
                self.ions_positions            = self.params.chain.get_positions()

        '''

        '''
        def set_potential(self, potential ):
                """ Take fz = V(x_0, y_0, z), the potential along Z axis 
                    and pass it to the Potential instance.

                """

                self.potential = potential

                if self.potential.config == 'generic':
                    fz             =  self.potential.fz_along_Z
                    self.potential.set_fz( fz )

                if not self.loaded:
                    self.load(self.chain)
        '''


        def set_used_electrodes(self, used_electrodes):

                self.potential.set_used_electrodes(used_electrodes)


        def set_electrode_voltages(self, electrode_voltages):
                
                self.potential.set_electrode_voltages(electrode_voltages)


        def set_trap_radial_freq(self, radial_freq, axis):
                """ Set radial_freq in Hertz

                """

                self.potential.set_trap_radial_freq(radial_freq, axis)


class Potential:

        def __init__(self, potential_config, trap_Y_radial_freq=1.8e6, trap_X_radial_freq=1.8e6,
                     optimization_parameters=OptimizationParameters(), 
                     max_fit_order = 3, z_expansion_order=4, X_fit_order=4, Y_fit_order=4):
                """ Contains information about total potential which is devided by 
                    0.5*M*omegaz**2 along the Z axis. 
                    func is a numpy.polynomial object, determining V(x_0,y_0, z), potential along the Z trap axis


                    trap_Y_radial_freq and trap_X_radial_freq are in Hertz
                """ 



                self.max_fit_order   = max_fit_order
                self.z_expansion_order = z_expansion_order
                self.Y_fit_order     = Y_fit_order
                self.X_fit_order     = X_fit_order

                self.set_axial_freq = NaN
                self.energy_deriv_along_Z   = NaN

                self.set_trap_radial_freq( trap_Y_radial_freq, 'Y')
                self.set_trap_radial_freq( trap_X_radial_freq, 'X')
                

                self.params          = optimization_parameters #simulation_parameters()

                self.config          = potential_config
                self.electrode_voltages_are_set = False



                if self.config == 'harmonic':
                        if axial_freq == 0:
                                raise Exception("Axial frequency is not set.")
                        self.axial_freq = self.params.axial_freq #Assuming that it is already multiplied by 2*np.pi before passing to Potential.


                #Have YOU not seen it through the end?                
                
        def set_fz(self, func):
                """ Set func which is a numpy.polynomial (poly1d) object, determining V(x_0,y_0, z), potential along the Z trap axis
                    Units of coefficients should be in S.I.

                """
                self.config = 'generic'

                self.energy_along_Z  =  self.params.charge *  func #We need to multiply by electron charge to get potential energy.
                

        def set_axial_freq(self, axial_freq ):

                self.axial_freq  = 2*np.pi * axial_freq


        def set_trap_radial_freq(self, radial_freq, axis):

                if axis == 'Y':
                    self.trap_Y_radial_freq  = 2*np.pi * radial_freq
                elif axis == 'X':
                    self.trap_X_radial_freq  = 2*np.pi * radial_freq



        def get_trap_radial_freq(self, axis):

                if axis == 'Y':
                    return self.trap_Y_radial_freq
                elif axis == 'X':
                    return self.trap_X_radial_freq 




        def set_used_electrodes(self, used_electrodes):

                self.used_electrodes = used_electrodes

        def set_electrode_voltages(self, electrode_voltages):

                if len(electrode_voltages) != len(self.used_electrodes):
                        self.electrode_voltages_are_set = False
                        raise Exception("Number of given electrode voltages does not match with the number of used_electrodes")

                else:

                    self.electrode_voltages = electrode_voltages
                    
                    self.pot            = self.get_potential()
                    self.V_along_Z_axis = self.get_V_along_Z_axis()
                    self.V_on_Y_axis_along_Z_axis = self.get_V_on_Y_axis_along_Z_axis()
                    self.V_on_X_axis_along_Z_axis = self.get_V_on_X_axis_along_Z_axis()
                    self.Y_1st_deriv, self.Y_2nd_deriv = self.get_Y_1st_and_2nd_deriv_fits_along_Z()
                    self.X_1st_deriv, self.X_2nd_deriv = self.get_X_1st_and_2nd_deriv_fits_along_Z()
                    self.fz_along_Z                    = self.get_fz_along_Z() 
                    print( 'fz', self.fz_along_Z )
                    self.set_fz(  self.fz_along_Z )

                    self.electrode_voltages_are_set = True


        #Brought in from Spectrum
        def get_potential(self):
                """ Return potential in 3 dimensions 

                """
                #Note: Substracting each electrode's potential at all points from its value at the trap center 
                #is good for fitting (so that even if fit with higher orders, the lower order coefficients of fit won't
                #change much).
                s = 0
                arr = 0
                for i in self.used_electrodes:
                    arr +=  self.electrode_voltages[s]*( self.params.pkl['EL_DC_{}'.format(i)] -  self.params.pkl['EL_DC_{}'.format(i)][self.params.trap_center] ) 
                    s += 1
                return arr 

                
        def get_V_along_Z_axis(self):
                """ Return an array of potentials along the Z axis with the given electrode voltages
                """
                return np.array( [self.pot[self.params.trap_center[0:2]+(z,)] for z in range(self.params.data_size[2])] )


        def get_V_on_Y_axis_along_Z_axis(self):
                """ Return a list of arrays of potentials on Y axis along the Z asis (one array
                for each Z position) with the given electrode voltages
                """
                arr = [ np.array( [self.pot[self.params.trap_center[0],y, z] for y in range(self.params.data_size[1])] ) for z in range(self.params.data_size[2]) ]
                return arr

        def get_V_on_X_axis_along_Z_axis(self):
                """ Return a list of arrays of potentials on X axis along the Z asis (one array
                for each Z position) with the given electrode voltages
                """
                return [ np.array( [ self.pot[x, self.params.trap_center[1], z] for x in range(self.params.data_size[0])] ) for z in range(self.params.data_size[2]) ]



        def get_Y_1st_and_2nd_deriv_fits_along_Z(self, Y_fit_order=4):
                """ Fit each array in V_on_Y_axis_along_Z_axis, and
                return an array of first and an array of second derivatives along 
                the Z axis with the given electrode voltages.

                """
                Y_1st_deriv  = []
                Y_2nd_deriv  = []
                
                for z in range(len(self.params.Z)):
                        fit_polynomial = polyfit( self.params.Y, self.V_on_Y_axis_along_Z_axis[z] , Y_fit_order)
                        
                        Y_1st_deriv.append( polyder(poly1d(fit_polynomial), 1)( self.params.Y[self.params.trap_center[1]] ) ) 
                        Y_2nd_deriv.append( polyder(poly1d(fit_polynomial), 2)( self.params.Y[self.params.trap_center[1]] ) ) 

                fit_Y_1st_deriv = poly1d( polyfit(  self.params.Z, Y_1st_deriv, self.z_expansion_order ) )
                fit_Y_2nd_deriv = poly1d( polyfit(  self.params.Z, Y_2nd_deriv, self.z_expansion_order ) )

                return fit_Y_1st_deriv, fit_Y_2nd_deriv


        def get_X_1st_and_2nd_deriv_fits_along_Z(self, X_fit_order=4):
                """ Fit each array in V_on_X_axis_along_Z_axis, and
                return an array of first and an array of second derivatives along 
                the Z axis with the given electrode voltages.

                """
                X_1st_deriv  = []
                X_2nd_deriv  = []
                
                for z in range(len(self.params.Z)):
                        fit_polynomial = polyfit( self.params.X, self.V_on_X_axis_along_Z_axis[z] , X_fit_order)
                        
                        X_1st_deriv.append( polyder(poly1d(fit_polynomial), 1)( self.params.X[self.params.trap_center[1]] ) ) 
                        X_2nd_deriv.append( polyder(poly1d(fit_polynomial), 2)( self.params.X[self.params.trap_center[1]] ) ) 

                fit_X_1st_deriv = poly1d( polyfit(  self.params.Z, X_1st_deriv, self.z_expansion_order ) )
                fit_X_2nd_deriv = poly1d( polyfit(  self.params.Z, X_2nd_deriv, self.z_expansion_order ) )

                return fit_X_1st_deriv, fit_X_2nd_deriv



        def get_fz_along_Z(self):
                """ Return f(x_0,y_0, z) in potential along Z axis.

                """

                return poly1d(polyfit( self.params.Z, self.V_along_Z_axis, self.z_expansion_order ))


        def get_Z_derivative(self, n):
                """ Return nth derivative of V(x_0,y_0, z) with respect to z (as a polynomial function of z).

                """

                return polyder(self.fz_along_Z, n)





        def Laplacian(self, zmin_index, zmax_index ):

                self.Z2 = self.get_Z_derivative(1)(self.potential.Z) 
                self.X2 = self.X_2nd_deriv(self.params.Z)
                self.Y2 = self.Y_2nd_deriv(self.params.Z)
                self.Laplace_eq = self.X2 + self.Y2 + self.Z2




        def get_radial_correction(self, radial):

                if radial == 'Y':
                        #Correction to Y radial frequency:
                        y2 = 0.5*poly1d(self.Y_2nd_deriv) #coeff of y2 in V(y,z) = y2*alpha*f(z) 
                        omega_y_correction = self.params.y_delta_coeff*y2 /(2*np.pi) #correction to the y radial frequency in Hertz, put z in mm
                        omega_y_correction_zgradient = 0.001*polyder(omega_y_correction,1) #in Hertz per microns unit, but you put z in mm

                        return omega_y_correction

                elif radial == 'X':
                        #Correction to X radial frequency:
                        x2 = 0.5*poly1d(spectrum.X_2nd_deriv) #coeff of y2 in V(y,z) = y2*alpha*f(z) 
                        omega_x_correction = self.params.x_delta_coeff*x2 /(2*np.pi) #correction to the x radial frequency in Hertz, put z in mm
                        omega_x_correction_zgradient = 0.001*polyder(omega_x_correction,1) #in Hertz per microns unit, but you put z in mm
           
                        return omega_x_correction


        def get_radial_frequency(self, radial):
                """ Get radial axis ('X' or 'Y') and return radial frequency along that direction after correction from
                static potential.

                """
                if radial == 'Y':
                        return self.trap_Y_radial_freq + self.get_radial_correction(radial) 
                elif radial == 'X':
                        return self.trap_X_radial_freq + self.get_radial_correction(radial) 



        def update(self, update_func):

                self.Vz += update_func
