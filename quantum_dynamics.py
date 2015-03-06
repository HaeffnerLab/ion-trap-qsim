#A collection of methods to find :
from __future__ import division, absolute_import, print_function, unicode_literals
import numpy as np 
import matplotlib.pylab as plt
import simulation_parameters
from numpy import absolute, real, linspace, arange, exp, pi
from scipy.misc import factorial  
from qutip import *
from operator_zoo import OperatorZoo

class Dynamics(OperatorZoo):
        
    def __init__(self, chain, lasers=[], pulses=[], time_precision=2.e-6):
 
        self.eta  =  0.05
        self.chain = chain
        self.time_precision = time_precision
        self.lasers = lasers
        for laser in self.lasers:
            if laser.ion_num > chain.num_of_ions or laser.ion_num < 1:
                raise Exception( "Laser ion number must be a number between 1 and {}.".format(chain.num_of_ions) )

        super(Dynamics, self).__init__()
        self.expectations = [[]]

        if not self.chain.motional_states_are_set:
            raise Exception("Motional states are not set.")
        self.chain.initialize_chain_electronic_states(lasers=lasers, pulses=lasers)

        self.chain_motion_hamiltonian = self.get_chain_motion_hamiltonian()
        self.construct_free_hamiltonian()
        #self.hamiltonian         =  self.get_chain_motion_hamiltonian() 
        
        
        #self.exponential         =  (-1.j*self.hamiltonian * self.time_precision ).expm() 
        #self.exponential         =  self.identity_op - 1.j*self.free_hamiltonian * self.time_precision 
        


    def simulate(self, times, example, ion, rf, DELTA, delta_radial_freq, omegax_gradient_at_100microns):

        if example == 5:
            psi0     = self.one_ion_excited_initial_state(ion)
            op_term1 = tensor(self.self_correlation_op_arr(ion))
            op = op_term1
            #modify local radial frequency of last two ions by self.delta_radial_freq:
            self.omegax[0][0] +=  self.delta_radial_freq
            self.omegax[self.N-1][self.N-1] += self.delta_radial_freq

            H = self.get_free_Hamiltonian() #+ self.AC_Stark_Hamiltonian(AC_Stark_ion, DELTA, OMEGA)                    
            
            # Find time evolution of psi0 and op
            output1 = mcsolve(H, psi0, times, [], [op])     


    def one_ion_excited_initial_state(self, ion):
        '''Return an initial state with all the ions in ground state 
        except one, given by ion, which is between 1 and self.N.
        '''
        arr = [0] * self.N
        arr[ion-1] = 1
        return self.ket(arr)



    def construct_propagator(self, times, time_independent=True, flag='all'):
        """Retrun propagator for Hamiltonian H. If flag=='all' return
        propagator at all time points given in times, if flag=='final' return
        propagator at time zero and times[-1].

        """
        if time_independent  == True:
            identity_op             =  (-1.j * self.hamiltonian * time_precision/100.).expm() 

            if flag == 'all':
                time_precision          =  times[1] - times[0]
                time_evol_op_list       =  []
                time_evol_op            =  identity_op
                exponential             =  (-1.j * self.hamiltonian * time_precision ).expm() 
                for t in enumerate(times):
                    time_evol_op_list      +=   [time_evol_op]
                    time_evol_op            =   exponential * time_evol_op

            elif flag == 'final':
                time_evol_op            =  identity_op #nntensor(qeye(2), qeye(2), qeye(M), qeye(M))
                exponential             =  (-1.j * self.hamiltonian * time_precision ).expm() 
                for t in enumerate(times):
                    time_evol_op        =   exponential * time_evol_op
                time_evol_op_list       =   [time_evol_op]
            
            self.time_evol_op_list  =  time_evol_op_list
            

    def construct_free_hamiltonian(self):

        self.free_hamiltonian  =  0

        #Add lasers Hamiltonian terms:
        for laser in self.lasers:
            self.free_hamiltonian += self.get_1st_order_hamiltonian(laser)

        #Add chain free Hamiltonian: 
        self.free_hamiltonian  +=  self.chain_motion_hamiltonian



    def get_laser_hamiltonian(self):

        H   =   0
        for laser in self.lasers:
            H += get_1st_order_hamiltonian(laser)

        return H


    def get_1st_order_hamiltonian(self, laser, regime='RWA'):
        """ Generate the first carrier(sideband) Hamiltonian in the Lamb-Dicke regime up to 
        1st order in eta**laser.sideband_num .

        """
        if regime == 'RWA':
            print("Simulation running in RWA regime")
            if laser.sideband_num > 0:
                op = ( 1.j * self.eta * exp(1.j*laser.phase) * self.a[laser.ion_num-1].dag() )**int(laser.sideband_num) / float(factorial(laser.sideband_num))
            elif laser.sideband_num < 0:     
                op = ( 1.j * self.eta * exp(1.j*laser.phase) * self.a[laser.ion_num-1] )**abs(int(laser.sideband_num)) / float(factorial(abs(laser.sideband_num)))
            else:
                op = 1

            op  =  self.sigma_plus[laser.ion_num-1] * op

            H   =  laser.intensity * ( op + op.dag() )  \
                    - self.chain.couplings[ laser.ion_num -1 ][ laser.ion_num -1 ] * sum( [ self.a[i].dag() * self.a[i] for i in range(self.chain.num_of_ions)  ] )



            #Add effect of laser detuning from local carrier or sideband frequency
            if laser.intensity != 0:
                H   +=  -1 * laser.detuning * self.sigmaz[laser.ion_num-1]

            return H


    def get_chain_motion_hamiltonian(self):
        """ Return the free evolution Hamiltonian of the chain motion in the regime where fast rotating terms 
        such as a^dag ** 2 are negligible.
        """

        s   =   0
        for i in range(self.chain.num_of_ions):
            for j in range(self.chain.num_of_ions):
                s  +=   self.chain.get_couplings()[i][j] * self.a[i].dag()*self.a[j]
        return s
    

    def get_first_blue_hamiltonian(self, ion_num, eta, omega, phase):
        return  1.j * eta*omega * ( exp(-1.j*phase) * self.sigma_plus[ion_num] * self.a[ion_num].dag() - exp(1.j*phase) * self.sigma_minus[ion_num] * self.a[ion_num] ) 


    def evolve_pure(self, time_interval, observables, time_dependent=False):
        psi0    =  self.chain.initial_state
        times = arange( time_interval[0], time_interval[1], self.time_precision )
        output1 = mcsolve(self.free_hamiltonian, psi0, times, [], observables)     

        self.expectations = output1.expect


    def evolve_free(self, time_interval, observables, time_dependent=False, state='pure'):

        self.expectations = [[]]
        self.time_evol_op = self.exponential

        try:
            self.chain.state         =  self.chain.initial_state
            times = arange( time_interval[0], time_interval[1], self.time_precision )
            if state == 'pure':
                for time in times:
                    self.expectations[0]  += [ real( expect( observables[0], self.time_evol_op * self.chain.state ) ) ]
                    self.chain.state       =  self.time_evol_op * self.chain.state
                 
            elif state == 'mixed':
                print("To be coded!")

        except AttributeError as e:
            print(e, "Chain initial state is not set.")


    #def do_measurement( observables ):

    #    for i, e in enumerate(observables):
    #        self.expectations  = real( expect( e, self.time_evol_op * rho0 * self.time_evol_op.dag() ) ) 
        


    def propagator_gen(self, time_precision):
        self.time_evol_op        =  self.identity_op
        exponential         =  (-1.j*self.hamiltonian * time_precision ).expm() 
        while True:
            yield self.time_evol_op
            self.time_evol_op        =   exponential * self.time_evol_op
            


    #def add_pulse( self, ion_num, pulse_duration, sideband_num = 1):
    #    self.pulse_hamiltonian = {}
    #    self.pulse_hamiltonian[str(ion_num)] = 
