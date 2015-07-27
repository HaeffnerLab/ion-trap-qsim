from __future__ import division, absolute_import, print_function, unicode_literals
import numpy as np
import pickle


class SimulationParameters(object):

        def __init__(self):
        
            
                #ion parameters
                self.mass = 40 * 1.6605402e-27 #40 amu in kg
                
                self.charge        = 1.6e-19
                self.e             = 1.6e-19
                
                self.coulomb_coeff = 2.30707955552e-28 # k =  U.e**2 / (4.0 * U.pi * U.eps0) 

                self.hbar = 1.05457266913e-34
                

                self.do_print = False

class PotentialData(object):


        def __init__(self):

                file_path          = '/home/trxw/Documents/dfr/codes/Fitting/BEMSolver_cfiles/NewCfiles/1microns_Center_-16_107_1270/data.pkl'
                self.pkl           = pickle.load(open(file_path, 'r'))
                self.data_size     = (31,31,331)
                self.trap_center   = (15, 15, 165)
                self.position_units = 1.e-3 #millimeter

                self.X = (self.pkl.X - self.pkl.X[self.trap_center[0]])*self.position_units
                self.Y = (self.pkl.Y - self.pkl.Y[self.trap_center[1]])*self.position_units
                self.Z = (self.pkl.Z - self.pkl.Z[self.trap_center[2]])*self.position_units
                
                

            

class OptimizationParameters(SimulationParameters, PotentialData):

        def __init__(self): #y_radial_freq_range=linspace(1.8e6, 2.e6, 10), x_radial_freq_range=linspace(1.8e6, 2.e6, 10)):


                SimulationParameters.__init__(self)
                PotentialData.__init__(self)

                #Correct the following lines:                
                
                
                #self.number_of_ions= number_of_ions

                #self.axial_freq    = 2*np.pi* 200.e3 #This value will be used in a Potential instance if potential_config is set to 'harmonic'

                