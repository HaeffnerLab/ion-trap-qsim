import pickle
from numpy import *
from scipy.optimize import basinhopping


import os,sys,inspect

import time
from scipy import linalg as LA 
from scipy import optimize
import numpy as np
from simulation_parameters import SimulationParameters, OptimizationParameters
from numpy.random import random
from ions import Chain
from ion_trap import *
import matplotlib.pylab as plt





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



class Optimization:


        def __init__(self, trap, spectral_seperation_ratio=.1, laser_orientation='Y', 
                     Y_field_max=0.02 * 1.e3 , X_field_max=0.02 * 1.e3, Z_field_max=0.2 * 1.e3,
                      voltage_scale=100., diff_voltage_scale=40., max_voltage_scale = 50.,
                      plot_cost = True, electrodes_connected=True):
            
            self.trap        = trap
            self.Y_field_max = Y_field_max

            self.X_field_max = X_field_max
            self.Z_field_max = Z_field_max
            
            self.voltage_scale = voltage_scale
            self.max_voltage_scale = max_voltage_scale
            self.diff_voltage_scale = diff_voltage_scale


            #self.initial_potential = initial_potential
            self.chain       =  trap.chain
            self.laser_orientation = laser_orientation
            self.spectral_seperation_ratio = spectral_seperation_ratio
            self.optimal_voltages_and_spectra = []
            self.cost = 1
            self.run = 1
            self.file = open('Optimization/Results/global_optimization_results_2modes_isolation.txt', 'w')

            self.electrodes_connected = electrodes_connected
            self.plot_cost = plot_cost
            self.plot_counter = 1

        def cost_function1(self, electrode_voltages):

            try:
                self.trap.set_electrode_voltages( electrode_voltages )
            except Exception:
                #self.randomize_voltages(electrode_voltages)
                self.cost += .1*abs(self.cost)
                return self.cost

            #Following condition must turn into a zigzag avoiding condition:
            if max(self.trap.chain.get_positions())-min(self.trap.chain.get_positions())< self.trap.chain.num_of_ions * 2.e-6:
                #self.randomize_voltages(electrode_voltages)
                self.cost += .4*abs(self.cost)

            else:

                normal_modes_with_radial_correction = self.trap.get_spectrum(self.laser_orientation)[1]
                mode_splitting_ratios = []
                for i in range(3):
                    mode_splitting_ratios.append( (normal_modes_with_radial_correction[i]-normal_modes_with_radial_correction[i+1])/(normal_modes_with_radial_correction[i+1] \
                        - normal_modes_with_radial_correction[i+2]) )
                            
                smallest_mode_splitting_ratio = min(mode_splitting_ratios)

                #Set cost criteria for spectrum:
                spectrum_cost = smallest_mode_splitting_ratio/self.spectral_seperation_ratio

                Efield_cost = (1./3) *( sum(  (  self.trap.potential.Y_1st_deriv( self.chain.get_positions() )/self.Y_field_max  )**2  )/self.chain.num_of_ions  +  \
                            sum(  (  self.trap.potential.X_1st_deriv( self.chain.get_positions() )/self.X_field_max  )**2  )/self.chain.num_of_ions  +\
                            sum(  (  self.trap.potential.Z_1st_deriv( self.chain.get_positions() )/self.Z_field_max  )**2  )/self.chain.num_of_ions  )                        


                voltage_cost = sum(  (                  self.trap.potential.electrode_voltages/self.voltage_scale  )**2  )/len( self.trap.potential.used_electrodes )

                self.cost = spectrum_cost + Efield_cost + voltage_cost

                print("\nCost is %.15f" % self.cost)
                print("Voltages are ", electrode_voltages)

            return self.cost



        def cost_function_new(self, electrode_voltages):

            try:
                self.trap.set_electrode_voltages( electrode_voltages )
            except Exception:
                #self.randomize_voltages(electrode_voltages)
                self.cost += 5.#.1*abs(self.cost)
                print "Dark region, cost: "+str(self.cost)
                return self.cost

            #Following condition must turn into a zigzag avoiding condition:
            if max(self.trap.chain.get_positions())-min(self.trap.chain.get_positions())< self.trap.chain.num_of_ions * 2.e-6:
                #self.randomize_voltages(electrode_voltages)
                self.cost += 5.#.1*abs(self.cost)
                print "Dark region, cost: "+str(self.cost)
                return self.cost

            else:

                normal_modes_with_radial_correction = self.trap.get_spectrum(self.laser_orientation)[1]
                
                '''
                mode_splitting_ratios = []
                for i in range(3):
                    mode_splitting_ratios.append( (normal_modes_with_radial_correction[i]-normal_modes_with_radial_correction[i+1])/(normal_modes_with_radial_correction[i+1] \
                        - normal_modes_with_radial_correction[i+2]) )
                            
                smallest_mode_splitting_ratio = min(mode_splitting_ratios)

                '''

                first_mode_splitting_ratio = (normal_modes_with_radial_correction[0]-normal_modes_with_radial_correction[1])/(normal_modes_with_radial_correction[1] \
                        - normal_modes_with_radial_correction[2]) 


                #Set cost criteria for spectrum:
                spectrum_cost = exp(first_mode_splitting_ratio/self.spectral_seperation_ratio)


                Efield_cost = (1./3) *( sum(  (  self.trap.potential.Y_1st_deriv( self.chain.get_positions() )/self.Y_field_max  )**2  )/self.chain.num_of_ions  +  \
                            sum(  (  self.trap.potential.X_1st_deriv( self.chain.get_positions() )/self.X_field_max  )**2  )/self.chain.num_of_ions  +\
                            5. *sum(  (  self.trap.potential.Z_1st_deriv( self.chain.get_positions() )/self.Z_field_max  )**2  )/self.chain.num_of_ions  )                        


                voltage_cost = sum(  (                  self.trap.potential.electrode_voltages/self.voltage_scale  )**2  )/len( self.trap.potential.used_electrodes )

                self.cost = spectrum_cost + Efield_cost + voltage_cost

                print("\nCost is %.15f" % self.cost)
                print"first_mode_splitting_ratio " + str(first_mode_splitting_ratio)

                print("Voltages are ", electrode_voltages)


                if self.plot_counter == 1 and self.plot_cost:
                    plt.ylim(0, 1.3*self.cost)
                    plt.ion()
                    plt.show()            
                self.plot_realtime( self.cost )

            return self.cost



        def cost_function_exp_voltages_punish(self, electrode_voltages):
            try:
                self.trap.set_electrode_voltages( electrode_voltages )
            except Exception:
                #self.randomize_voltages(electrode_voltages)
                self.cost += .01*abs(self.cost)
                return self.cost

            #Following condition must turn into a zigzag avoiding condition:
            if max(self.trap.chain.get_positions())-min(self.trap.chain.get_positions())< self.trap.chain.num_of_ions * 2.e-6:
                #self.randomize_voltages(electrode_voltages)
                self.cost += .01*abs(self.cost)
                return self.cost
            else:

                normal_modes_with_radial_correction = self.trap.get_spectrum(self.laser_orientation)[1]
                
                '''
                mode_splitting_ratios = []
                for i in range(3):
                    mode_splitting_ratios.append( (normal_modes_with_radial_correction[i]-normal_modes_with_radial_correction[i+1])/(normal_modes_with_radial_correction[i+1] \
                        - normal_modes_with_radial_correction[i+2]) )
                            
                smallest_mode_splitting_ratio = min(mode_splitting_ratios)

                '''

                first_mode_splitting_ratio = (normal_modes_with_radial_correction[0]-normal_modes_with_radial_correction[1])/(normal_modes_with_radial_correction[1] \
                        - normal_modes_with_radial_correction[2]) 


                #Set cost criteria for spectrum:
                spectrum_cost = first_mode_splitting_ratio/self.spectral_seperation_ratio


                Efield_cost = (1./3) *( sum(  (  self.trap.potential.Y_1st_deriv( self.chain.get_positions() )/self.Y_field_max  )**2  )/self.chain.num_of_ions  +  \
                            sum(  (  self.trap.potential.X_1st_deriv( self.chain.get_positions() )/self.X_field_max  )**2  )/self.chain.num_of_ions  +\
                            sum(  (  self.trap.potential.Z_1st_deriv( self.chain.get_positions() )/self.Z_field_max  )**2  )/self.chain.num_of_ions  )                        


                voltage_cost = sum(     exp( (self.trap.potential.electrode_voltages/self.voltage_scale  )**2 )   )/len( self.trap.potential.used_electrodes )

                self.cost = spectrum_cost + Efield_cost + voltage_cost


                print "\nCost is %.15f" % self.cost
                print "Voltages: " + str( list(electrode_voltages) )


                #self.save( self.cost, electrode_voltages )
                self.run += 1

            return self.cost



        def cost_function_voltages_and_max_voltages_exp_punish(self, electrode_voltages):


            try:
                self.trap.set_electrode_voltages( electrode_voltages )

            except Exception:
                #self.randomize_voltages(electrode_voltages)
                self.cost += .001*abs(self.cost)
                print "\nBad Cost is %.15f" % self.cost
                return self.cost
                #raise Exception("Whatever")

            #Following condition must turn into a zigzag avoiding condition:
            if max(self.trap.chain.get_positions())-min(self.trap.chain.get_positions())< self.trap.chain.num_of_ions * 1.e-6:
                #self.randomize_voltages(electrode_voltages)
                self.cost += .001*abs(self.cost)
                print "\nChain crushed. Cost is %.15f" % self.cost
                return self.cost
                #raise Exception("\nChain Crushed")


            else:











                normal_modes_with_radial_correction = self.trap.get_spectrum(self.laser_orientation)[1]
                
                '''
                mode_splitting_ratios = []
                for i in range(3):
                    mode_splitting_ratios.append( (normal_modes_with_radial_correction[i]-normal_modes_with_radial_correction[i+1])/(normal_modes_with_radial_correction[i+1] \
                        - normal_modes_with_radial_correction[i+2]) )
                            
                smallest_mode_splitting_ratio = min(mode_splitting_ratios)

                '''

                first_mode_splitting_ratio = (normal_modes_with_radial_correction[0]-normal_modes_with_radial_correction[1])/(normal_modes_with_radial_correction[1] \
                        - normal_modes_with_radial_correction[2]) 


                #Set cost criteria for spectrum:
                spectrum_cost = 0*first_mode_splitting_ratio/self.spectral_seperation_ratio
                spectrum_exp_cost = 10*exp( (first_mode_splitting_ratio/0.5)**2 ) -1

                Efield_cost = (1./3) *( sum(  (  self.trap.potential.Y_1st_deriv( self.chain.get_positions() )/self.Y_field_max  )**2  )/self.chain.num_of_ions  +  \
                            sum(  (  self.trap.potential.X_1st_deriv( self.chain.get_positions() )/self.X_field_max  )**2  )/self.chain.num_of_ions  +\
                            0*sum(  exp( (  self.trap.potential.Z_1st_deriv( self.chain.get_positions() )/self.Z_field_max  )**2)  )/self.chain.num_of_ions  )                        


                voltage_cost = 0* sum(     (self.trap.potential.electrode_voltages/self.voltage_scale  )**2    )/len( self.trap.potential.used_electrodes )
                #diff_voltages_cost = sum(     exp( (diff(electrode_voltages)/self.diff_voltage_scale  )**2 )   )/len( self.trap.potential.used_electrodes - 1 )
               
                #Punish max voltage significantly if it was above 80
                max_voltage_cost = 0 #exp( (max(abs(self.trap.potential.electrode_voltages))/self.max_voltage_scale)**2 )

                #max_voltages_cost  = (max(abs(self.trap.potential.electrode_voltages)) - min(abs(self.trap.potential.electrode_voltages)))/self.diff_voltage_scale


                self.cost = spectrum_cost + spectrum_exp_cost + 0*Efield_cost + voltage_cost + max_voltage_cost


                print "\nOK Cost is %.15f" % self.cost
                print "Voltages: " + str( list(electrode_voltages) )


                #self.save( self.cost, electrode_voltages )
                print "\nfirst_mode_splitting_ratio: %.5f" % first_mode_splitting_ratio
                
                last_ok_cost = self.cost

                if self.plot_counter == 1 and self.plot_cost:
                    plt.ylim(0, 1.3*self.cost)
                    plt.ion()
                    plt.show()            
                self.plot_realtime( self.cost )
                self.run += 1
                
                return self.cost




        def cost_function_play(self, electrode_voltages):

            
            self.trap.set_electrode_voltages( electrode_voltages )

            #Following condition must turn into a zigzag avoiding condition:
            if max(self.trap.chain.get_positions())-min(self.trap.chain.get_positions())< self.trap.chain.num_of_ions * 2.e-6:
                #self.randomize_voltages(electrode_voltages)
                self.cost = 1000.#*abs(self.cost)
                print "\nChain crushed. Cost is %.15f" % self.cost
                #raise Exception("\nChain Crushed")


            else:


                normal_modes_with_radial_correction = self.trap.get_spectrum(self.laser_orientation)[1]
                
                '''
                mode_splitting_ratios = []
                for i in range(3):
                    mode_splitting_ratios.append( (normal_modes_with_radial_correction[i]-normal_modes_with_radial_correction[i+1])/(normal_modes_with_radial_correction[i+1] \
                        - normal_modes_with_radial_correction[i+2]) )
                            
                smallest_mode_splitting_ratio = min(mode_splitting_ratios)

                '''

                first_mode_splitting_ratio = (normal_modes_with_radial_correction[0]-normal_modes_with_radial_correction[1])/(normal_modes_with_radial_correction[1] \
                        - normal_modes_with_radial_correction[2]) 


                #Set cost criteria for spectrum:
                spectrum_cost = first_mode_splitting_ratio/self.spectral_seperation_ratio
                spectrum_exp_cost = exp( (first_mode_splitting_ratio/0.1)**2 ) -1

                Efield_cost = (1./3) *( sum(  (  self.trap.potential.Y_1st_deriv( self.chain.get_positions() )/self.Y_field_max  )**2  )/self.chain.num_of_ions  +  \
                            sum(  (  self.trap.potential.X_1st_deriv( self.chain.get_positions() )/self.X_field_max  )**2  )/self.chain.num_of_ions  +\
                            sum(  exp( (  self.trap.potential.Z_1st_deriv( self.chain.get_positions() )/self.Z_field_max  )**2)  )/self.chain.num_of_ions  )                        


                voltage_cost = 0* sum(     (self.trap.potential.electrode_voltages/self.voltage_scale  )**2    )/len( self.trap.potential.used_electrodes )
                #diff_voltages_cost = sum(     exp( (diff(electrode_voltages)/self.diff_voltage_scale  )**2 )   )/len( self.trap.potential.used_electrodes - 1 )
               
                #Punish max voltage significantly if it was above 80
                max_voltage_cost = 0 #exp( (max(abs(self.trap.potential.electrode_voltages))/self.max_voltage_scale)**2 )

                #max_voltages_cost  = (max(abs(self.trap.potential.electrode_voltages)) - min(abs(self.trap.potential.electrode_voltages)))/self.diff_voltage_scale


                self.cost = spectrum_cost + spectrum_exp_cost + Efield_cost + voltage_cost + max_voltage_cost


                print "\nCost is %.15f" % self.cost
                print "Voltages: " + str( list(electrode_voltages) )


                #self.save( self.cost, electrode_voltages )
                print "\nfirst_mode_splitting_ratio: %.5f" % first_mode_splitting_ratio
                
                last_ok_cost = self.cost

                if self.plot_counter == 1 and self.plot_cost:
                    plt.ylim(0, 1.3*self.cost)
                    plt.ion()
                    plt.show()            
                self.plot_realtime( self.cost )
                self.run += 1
                
                return self.cost










        def save(self, cost, electrode_voltages):
            self.file.write( "\n############################################################" )
            self.file.write( "\nRun %i " % self.run )
            self.file.write( "\nCost %.10f" % cost )
            self.file.write( "\nVoltages: \n" + str( list(electrode_voltages) ) )



        def plot_realtime(self, point):



            plt.scatter(self.plot_counter, point)
            plt.draw()
            self.plot_counter += 1




        

        def randomize( self, electrode_voltages):
            pass