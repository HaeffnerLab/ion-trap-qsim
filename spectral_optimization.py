import pickle
from numpy import *
from scipy.optimize import basinhopping


import os,sys,inspect


from scipy import linalg as LA 
from scipy import optimize
import numpy as np
from simulation_parameters import SimulationParameters
from numpy.random import random
from ions import Chain
from ion_trap import *
import matplotlib.pylab as plt
import seaborn as sns
sns.set_context('poster')

#%matplotlib inline 


class Spectrum:

        def __init__(self, optimization_parameters, potential, initial_positions_guess=[]):
                #used_electrodes and electrode_voltages should be in the same order.

                self.params             =  optimization_parameters
                
                self.potential          =   potential                
                

                #self.set_ions_positions()
                #if there was a solution to minimum potential, then:
                #self.get_spectrum()



        def does_trap(self):
                """Check to see if the z dependence of potential allows trapping
                """
                pass



        def get_equilibrium_positions(self):
                
                self.params.dummy_trap.get_positions(number_of_ions, potential, self.initial_positions_guess)

        

        def get_spectrum(self, laser_orientation='Y'):

                
                couplings             = self.params.chain.get_couplings()

                #Having a 729 laser with general angle could be added later:
                couplings_with_radial_correction = couplings
                if laser_orientation   == 'Y':
                    couplings_with_radial_correction  +=  self.get_radial_correction('Y')(self.ions_positions)
                elif laser_orientation == 'X':
                    couplings_with_radial_correction  +=  self.get_radial_correction('X')(self.ions_positions)



                normal_modes_from_Vz  = np.sort(abs(LA.eig(couplings)[0]))[::-1]  #This is the normal modes due to positioning of ions
                normal_modes_with_radial_correction = np.sort(abs(LA.eig(couplings_with_radial_correction)[0]))[::-1]  #This is the normal modes due to a generic potential

                return normal_modes_from_Vz, normal_modes_with_radial_correction





class SimulatedAnnealing:
        def __init__(self, initial_electrode_voltages, number_of_ions, spectral_seperation_ratio=0.01, voltage_scale=40.):

            self.number_of_ions   =  number_of_ions 
            self.initial_electrode_voltages   = initial_electrode_voltages
            self.optimization_parameters    = OptimizationParameters(number_of_ions, y_radial_freq,
                                                                x_radial_freq)

            self.spectral_seperation_ratio  = spectral_seperation_ratio
            self.voltage_scale              = voltage_scale

            sp   =  Spectrum(self.optimization_parameters, self.initial_electrode_voltages, initial_positions_guess)
            self.initial_positions_guess = []
            #self.sim_temperature = sim_temperature

            self.optimal_voltages_and_spectra = []



        def displacement_probability(self, voltages):
            pass

        def run(self, num_of_runs):

            init_voltages_guess = self.initial_electrode_voltages
            for i in range(num_of_runs):
                    displaced_voltages = next_point( init_voltages, 0.05  )
                    
                    sp = Spectrum(self.optimization_parameters, displaced_voltages, initial_positions_guess)
                    if cost_is_lower( init_voltages + displacement,  ): 
                            init_voltages += displacement
            
            position_list.append(init)



        def next_point(init_vec, mesh):
        
            init_vec += (2*random(len(init_vec))-1.)*mesh

            return init_vec

    
        def cost_is_lower(self, spectrum, critical_ratio=.1):
            

            try:
                    spectrum.set_ions_positions()
            except ValueError, e:
                    return False

            normal_modes_with_radial_correction = spectrum.get_spectrum(electrode_voltages)[1]
            for i in range(3):
                if (normal_modes_with_radial_correction[i]-normal_modes_with_radial_correction[i+1])/(normal_modes_with_radial_correction[i+1] \
                    - normal_modes_with_radial_correction[i+2]) <= critical_ratio:
                    
                    save(electrode_voltages, spectrum)
                    return True
            
            return False


        def cost_function(self, normal_modes_with_radial_correction, potential ):
            r = ( (normal_modes_with_radial_correction[0]-normal_modes_with_radial_correction[1])/(normal_modes_with_radial_correction[1] \
                    -normal_modes_with_radial_correction[2]) )/self.spectral_seperation_ratio

            return r +  sum(potential.Y_1st_deriv( ions_positions )/self.Y_field_max)/self.ions_positions  +  \
                        sum(potential.X_1st_deriv( ions_positions )/self.X_field_max)/self.ions_positions  +  \
                        sum((potential.electrode_voltages/self.voltage_scale)**2)


        def save(self, electrode_voltages, sp):

            self.optimal_voltages_and_spectra.append([electrode_voltages, sp])




