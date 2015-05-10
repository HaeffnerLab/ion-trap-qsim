import pickle
from numpy import *
import numpy as np


import os,sys,inspect


from scipy import linalg as LA 
from scipy import optimize
import numpy as np
from simulation_parameters import simulation_parameters
p = simulation_parameters()
from ions import Chain
from ion_trap import IonTrap
import matplotlib.pylab as plt
import seaborn as sns
sns.set_context('poster')

#%matplotlib inline 

class Parameters:
        def __init__(self):
                self.mass = 6.64e-26
                self.charge = 1.6e-19
                self.y_radial_freq = 2*np.pi * 1.9e6
                self.x_radial_freq = 2*np.pi * 1.9e6
                self.y_delta_coeff = 1.e6*self.charge/(2*self.y_radial_freq*self.mass) #conversion to frequency and to turn y2 derivative in SI units
                self.x_delta_coeff = 1.e6*self.charge/(2*self.x_radial_freq*self.mass) #conversion to frequency and to turn y2 derivative in SI units
                self.omegax = 2.0e6
                self.omegay = 2.0e6
                self.init_omegaz = 100.e3
                self.N = 15
                self.chain = Chain(N, 2)
                self.dummy_trap = IonTrap( omegax , omegaz)
self.dummy_trap.load(chain)
self.init_zpositions =   np.array( chain.get_positions() )


class Spectrum:

        def __init__(self, electrode_voltages, used_electrodes, pickle_file, data_size=(21,21,141), trap_center=(10, 10, 75), max_fit_order=3, z_expansion_order=4, Y_fit_order=4, X_fit_order=4):
                #used_electrodes and electrode_voltages should be in the same order.


                self.max_fit_order = max_fit_order
                self.z_expansion_order = z_expansion_order
                self.electrode_voltages = electrode_voltages
                self.used_electrodes = used_electrodes
                self.trap_center = trap_center

                
                self.pkl = pickle_file
                #print self.trap_center
                self.X = self.pkl.X - self.pkl.X[self.trap_center[0]]
                self.Y = self.pkl.Y - self.pkl.Y[self.trap_center[1]]
                self.Z = self.pkl.Z - self.pkl.Z[self.trap_center[2]]
                self.data_size = data_size
                
                self.pot = self.get_potential()
                self.V_along_Z_axis = self.get_V_along_Z_axis()
                self.V_on_Y_axis_along_Z_axis = self.get_V_on_Y_axis_along_Z_axis()
                self.V_on_X_axis_along_Z_axis = self.get_V_on_X_axis_along_Z_axis()
                self.Y_1st_deriv, self.Y_2nd_deriv = self.get_Y_1st_and_2nd_deriv_fits_along_Z()
                self.X_1st_deriv, self.X_2nd_deriv = self.get_X_1st_and_2nd_deriv_fits_along_Z()

                self.fz_along_Z               = self.get_fz_along_Z() 
                self.params   =  Parameters()

        def get_potential(self):
                """ Return potential in 3 dimensions 

                """
                #Note: Substracting each electrode's potential at all points from its value at the trap center 
                #is good for fitting (so that even if fit with higher orders, the lower order coefficients of fit won't
                #change much).
                s = 0
                arr = 0
                for i in self.used_electrodes:
                    arr +=  self.electrode_voltages[s]*( self.pkl['EL_DC_{}'.format(i)] -  self.pkl['EL_DC_{}'.format(i)][self.trap_center] ) 
                    s += 1
                return arr 
                
        def get_V_along_Z_axis(self):
                """ Return an array of potentials along the Z axis with the given electrode voltages
                """
                return np.array( [self.pot[self.trap_center[0:2]+(z,)] for z in range(self.data_size[2])] )


        def get_V_on_Y_axis_along_Z_axis(self):
                """ Return a list of arrays of potentials on Y axis along the Z asis (one array
                for each Z position) with the given electrode voltages
                """
                arr = [ np.array( [self.pot[self.trap_center[0],y, z] for y in range(self.data_size[1])] ) for z in range(self.data_size[2]) ]
                return arr

        def get_V_on_X_axis_along_Z_axis(self):
                """ Return a list of arrays of potentials on X axis along the Z asis (one array
                for each Z position) with the given electrode voltages
                """
                return [ np.array( [ self.pot[x, self.trap_center[1], z] for x in range(self.data_size[0])] ) for z in range(self.data_size[2]) ]



        def get_Y_1st_and_2nd_deriv_fits_along_Z(self, Y_fit_order=4):
                """ Fit each array in V_on_Y_axis_along_Z_axis, and
                return an array of first and an array of second derivatives along 
                the Z axis with the given electrode voltages.

                """
                Y_1st_deriv  = []
                Y_2nd_deriv  = []
                
                for z in range(len(self.Z)):
                        fit_polynomial = polyfit( self.Y, self.V_on_Y_axis_along_Z_axis[z] , Y_fit_order)
                        
                        Y_1st_deriv.append( polyder(poly1d(fit_polynomial), 1)( self.Y[self.trap_center[1]] ) ) 
                        Y_2nd_deriv.append( polyder(poly1d(fit_polynomial), 2)( self.Y[self.trap_center[1]] ) ) 

                fit_Y_1st_deriv = poly1d( polyfit(  self.Z, Y_1st_deriv, self.z_expansion_order ) )
                fit_Y_2nd_deriv = poly1d( polyfit(  self.Z, Y_2nd_deriv, self.z_expansion_order ) )

                return fit_Y_1st_deriv, fit_Y_2nd_deriv


        def get_X_1st_and_2nd_deriv_fits_along_Z(self, X_fit_order=4):
                """ Fit each array in V_on_X_axis_along_Z_axis, and
                return an array of first and an array of second derivatives along 
                the Z axis with the given electrode voltages.

                """
                X_1st_deriv  = []
                X_2nd_deriv  = []
                
                for z in range(len(self.Z)):
                        fit_polynomial = polyfit( self.X, self.V_on_X_axis_along_Z_axis[z] , X_fit_order)
                        
                        X_1st_deriv.append( polyder(poly1d(fit_polynomial), 1)( self.X[self.trap_center[1]] ) ) 
                        X_2nd_deriv.append( polyder(poly1d(fit_polynomial), 2)( self.X[self.trap_center[1]] ) ) 

                fit_X_1st_deriv = poly1d( polyfit(  self.Z, X_1st_deriv, self.z_expansion_order ) )
                fit_X_2nd_deriv = poly1d( polyfit(  self.Z, X_2nd_deriv, self.z_expansion_order ) )

                return fit_X_1st_deriv, fit_X_2nd_deriv



        def get_fz_along_Z(self):
                """ Return f(x_0,y_0, z) in potential along Z axis.

                """

                return poly1d(polyfit( self.Z, self.V_along_Z_axis, self.z_expansion_order ))


        def get_Z_derivative(self, n):
                """ Return nth derivative of V(x_0,y_0, z) with respect to z (as a polynomial function of z).

                """

                return polyder(self.fz_along_Z, n)



        def does_trap(self):
                """Check to see if the z dependence of potential allows trapping
                """
                self.fz_along_Z


        def get_equilibrium_positions(self):
                if does_trap():
                        pass

        

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


        def get_spectrum(self):

                pass

'''
#harmonic potential couplings:
couplings_harmonic = chain.get_couplings()#
print couplings_harmonic.shape
#Correction for x^2 * z^2 term:
omegax_gradient_at_100microns = +2*np.pi * (1)*1.e3/(1.e-6)
def get_radial_freq_corrections(omegax_gradient_at_100microns):
    
    delta_omegax = 0
    
    delta_omegax += 2*np.pi*omega_correction#omegax_gradient_at_100microns * zpositions**2/(2*100.e-6)

    return delta_omegax 


couplings_with_x2z2_correction = np.array( [ row for row in couplings_harmonic ])

radial_corrections = get_radial_freq_corrections(omegax_gradient_at_100microns)

for i in range(N):
    couplings_with_x2z2_correction[i][i] += radial_corrections[i]

normal_modes_harmonic = np.sort(abs(LA.eig(chain.get_couplings())[0])/(2*np.pi*1.e6))[::-1]
normal_modes_with_x2z2_correction = np.sort(abs(LA.eig(couplings_with_x2z2_correction)[0])/(2*np.pi*1.e6))[::-1]
                                    
plt.plot(normal_modes_harmonic, 'rs', normal_modes_with_x2z2_correction , 'bs')
plt.show()
'''





class inspector:

        def __init__(self, spectrum):
            self.spectrum  =  spectrum
            self.Z2 = self.spectrum.get_Z_derivative(1)(self.spectrum.Z) 
            self.X2 = self.spectrum.X_2nd_deriv(self.spectrum.Z)
            self.Y2 = self.spectrum.Y_2nd_deriv(self.spectrum.Z)

        def Laplace_eq(self):

            print self.X2 + self.Y2 + self.Z2





class optimize:
        def __init__(self, init_voltages):
                self.init_guess = init_voltages











'''





def get_Y2fz(self):


def get_X2fz(self):


def get_Y_field(electrode_voltages):
        """ Return an array of E_y electric fields along the Z asis with the given electrode voltages
        """
'''


if __name__ == '__main__':

        file_path = '/home/trxw/Documents/dfr/codes/Fitting/BEMSolver_cfiles/NewCfiles/Cfiles_ZCenter1240_copy/newdata.pkl'
        pkl         = pickle.load(open(file_path, 'r'))
        
        used_electrodes = range(1,24)
        
        '''
        electrode_voltages = array([  0.        ,  13.37686788,  18.38449342, -22.15026567,
          2.93369668,  -5.86005482,   6.45600753, -37.56846609,
          9.8046246 ,   8.82554234,   7.44922786,   4.57791903,
         15.3129863 ,  26.69914406,  31.66854976, -17.14898273,
         -3.48055587,   2.89242092,   8.85400634,   9.52939129,
          9.20496382,   7.7399846 ,   4.28803344,  -0.18530892])
        
        

        electrode_voltages = array([  0.        ,   0.23600562,   0.53091654,  -0.74890046,
         -5.36853701,  18.79960142, -11.27842687,   2.08858257,
          0.97143379,   0.29995416,   0.1768433 ,   0.08985979,
          1.06790129,   3.01670797,   8.46884727,   1.25495873,
         -6.96273258,   3.40082774,   7.20873915,   2.39018   ,
          0.92874822,   0.50319264,   0.1901779 ,  -0.15647462])

        '''
        electrode_voltages = array([  0.,   0.,   0.,   0.,   5.,  0.,   0.,   0.,   0.,   0.,   0.,
         0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
         0.])
        
    


        spectrum = Spectrum(electrode_voltages, used_electrodes, pkl)
 

        #Ispect the result:
        #Is Laplace eq satisfied at all points along Z axis? 
        insp = inspector(spectrum)
        insp.Laplace_eq()

