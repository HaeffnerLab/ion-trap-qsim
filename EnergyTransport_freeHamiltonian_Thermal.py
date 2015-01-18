#!/usr/bin/env python

import matplotlib.pylab as plt
from qutip import *
from scipy import *
import numpy as np
from equilbrium_positions import equilibrium_positions
import simulation_parameters
from scipy import linalg as LA 


class chain():
	"""
	Create a chain with the given parameters
	"""

	def __init__(self, omegax, N, M, omega_x, omega_z, equal_distances):
		self.omegax = omegax
		self.omega_x = omega_x
		self.omega_z = omega_z
		self.N = N 
   		self.M = M
		self.equal_distances = equal_distances
    		# Create an element of local modes basis, arr is the array of local modes, 
    		# e.g., arr = [0, 2, 1] gives |0,2,1>: Note len(arr) must be = N and 
    		# max(arr) must be < M:
    		self.ket = lambda arr: tensor([basis(self.M, arr[i]) for i in range(len(arr))])
    		self.ann_op = lambda i: tensor( [ destroy(self.M) if j == i else qeye(self.M) for j in range(self.N) ] )
            
		self.local_n_fock = lambda i: tensor([ create(self.M)*destroy(self.M) if j == i else qeye(self.M) for j in range(self.N) ])	
		# Create local Fock space destruction operators:
	        self.a = [self.ann_op(i) for i in range(self.N)]
	
		self.fock_identity = tensor( [ qeye(self.M) for i in range(self.N) ] )
	
	def simulate(self, ions, times, example, ions_init_fock_state, ion = 1, DELTA=0, OMEGA=0):
	        N = self.N
        	M = self.M
        	omegax = self.omegax
                middle_states = ions_init_fock_state[1:self.N-1]             
		if example == 1: #initial state (|01>+|10>)/sqrt(2), op = |01><10| + |10><01|, for ions[0] and ions[1]   
			psi0 = self.DFS_initial_state(ions)
        		op_term1 = tensor(self.op_arr(ions))
			op = op_term1 + op_term1.dag()
			# Find time evolution of psi0 and op with free Hamiltonian
    			H = self.free_Hamiltonian()
    			output1 = mcsolve(H, psi0, times, [], [op])
		
                
    		elif example == 2: # This is just a special case of example 4
      			psi0 = self.one_ion_excited_initial_state(ion)
    			op_term1 = tensor(self.self_correlation_op_arr(ion))
    			op = op_term1
			# Find time evolution of psi0 and op with free Hamiltonian
    			H = self.free_Hamiltonian()
    			output1 = mcsolve(H, psi0, times, [], [op])
                
                elif example == 3:
        		psi0 = self.one_ion_entangled_state(ion)
        		op_term1 = tensor(self.one_ion_1to0_projector_op_arr(ion))
        		op = op_term1
-			# Find time evolution of psi0 and op with free Hamiltonian
    			H = self.free_Hamiltonian()
    			output1 = mcsolve(H, psi0, times, [], [op])


        	elif example == 4:
                        ion_init_fock_state = ions_init_fock_state[ion-1]
			psi0 = self.one_ion_excited_initial_state(ion_init_fock_state, ion)
    			op_term1 = sum( [(float(n)/ion_init_fock_state) * tensor(self.self_correlation_op_arr(n, ion)) for n in range(1, sum(ions_init_fock_state)+1)] )
    			op = op_term1
			H = self.free_Hamiltonian() 
	
			# Find time evolution of psi0 and op
    			output1 = mcsolve(H, psi0, times, [], [op])
			
		elif example == 5:
			op_term1 = tensor(self.op_arr(ions))
  			op = op_term1 + op_term1.dag()
			
			psi0 = self.DFS_initial_state(ions, middle_states)
        		# Find time evolution of psi0 and op with free Hamiltonian
    			H = self.free_Hamiltonian()
    			output1 = mcsolve(H, psi0, times, [], [op])
				 
	
        	elif example == 6: #Far detuned 2nd blue sideband
                        eta = .05
			sigma_plus = create(2) 
                        expo_ion1 = 1.j*eta * ( self.a[0] + (self.a[0]).dag())
                        expo_ion1_series = expo_ion1**2/2 + expo_ion1**3/6
			H0 = self.free_Hamiltonian() 
			Ht_term1 = 0.5*OMEGA * tensor(sigma_plus, expo_ion1_series) 
			Ht_term2 = Ht_term1.dag()
			H = [H0, [H1, 'exp(-1.j*DELTA*t)']]
			output1 = mcsolve(H, psi0, times, [], [op])
		
				
    			
    		#Plot and save simulation data:    
    		if example == 1: 
			plt.figure(example*self.N)
			plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
	        	plt.xlabel('Time ($\mu$s)')
			plt.ylabel('$C_{a}$')
	        	plt.ylim(-.5,1)
			#plt.title('Expectation value of O = |01><10|+|10><01| for ions {} in a chain of {} ions, with max = {}'.format(ions,					N, max(output1.expect[0])))
			plt.xlim(0, times[-1]*1.e6)
			#Save the graph:
			plt.savefig('{}ions_example{}_DFSions{}_EqualDistanes{}_Radial{:.2f}MHz_Axial{:.2f}KHz_TwoOuterIonsNeighborTunneling{:.2f}KHz.pdf'.format(self.N, 
                example, ions, str(self.equal_distances), self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
			#Save the data:
			qsave(output1, '{}ions_DFSions{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_TwoOuterIonsNeighborTunneling{:.2f}KHz_data'.format(self.N, ions, example,  str(self.equal_distances), self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )
		
		elif example == 2: 
			plt.figure(example*(self.N+25)) # 25 is the max number of runs
			plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
	        	plt.xlabel('Time ($\mu$s)')
			plt.ylabel('Probability')
 	        	plt.ylim(0,1)
			plt.title('Probability of finding ion number {} in excited state, in a chain of {} ions, with max = {}'.format(ion, 
				N, max(output1.expect[0])))
        		plt.xlim(0, times[-1]*1.e6)
			#Save the graph:
			plt.savefig('{}ions_ion{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg'.format(self.N, ion, example,  str(self.equal_distances),
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
			qsave(output1, '{}ions_ion{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz_data'.format(self.N, ion, example,  str(self.equal_distances),
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )

    		elif example == 3:
			plt.figure(example*(self.N+25))
			plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
	        	plt.xlabel('Time ($\mu$s)')
			plt.ylabel('Re Transition amplitude')
        		plt.ylim(-.5,.5)
  	        	plt.title('Transition amplitude between |0> and |1> motional states of ion number {}, initially in (|0>+|1>)/sqrt(2) state, in a chain of {} ions, with max = {}'.format(ion, N, max(absolute(output1.expect[0]))))
			plt.xlim(0, times[-1]*1.e6)
			#Save the graph:
			plt.savefig('{}ions_ion{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg'.format(self.N, ion, example, str(self.equal_distances),
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
			#Save the data:
			qsave(output1, '{}ions_ion{}_example{}_Radial{:.2f}MHz_EqualDistances{}_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz_data'.format(self.N, ion, example, self.omega_x/(2*np.pi* 1.e6), str(self.equal_distances),
					 self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )

		elif example == 4:
			plt.figure(example*(self.N+25)) # 25 is the max number of runs
			plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
	        	plt.xlabel('Time ($\mu$s)')
			plt.ylabel('<E>')
 	        	plt.ylim(0,1)
        		plt.xlim(0, times[-1]*1.e6)
			#Save the graph:
			plt.savefig('{}ions_EnergyTransport_Ion{}_IonsInitState{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_Ion1-2Tunneling{:.2f}KHz.pdf'.format(self.N, ion, ions_init_fock_state, str(self.equal_distances),
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
			qsave(output1, '{}ions_EnergyTransport_Ion{}_IonsInitState{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_Ion1-2Tunneling{:.2f}KHz.pdf'.format(self.N, ion, ions_init_fock_state, str(self.equal_distances),
                                         self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )
	

      
		elif example == 5:
                        plt.figure(example*self.N)
                        plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
                        plt.xlabel('Time ($\mu$s)')
                        plt.ylabel('$C_{a}$')
                        plt.ylim(-1,1)
                        #plt.title('Expectation value of O = |01><10|+|10><01| for ions {} in a chain of {} ions, with max = {}'.format(ions,                                   N, max(output1.expect[0])))
                        plt.xlim(0, times[-1]*1.e6)
                        #Save the graph:
                        plt.savefig('{}ions_example{}_DFSions{}_MiddIonsInitState{}_EqualDistanes{}_Radial{:.2f}MHz_Axial{:.2f}KHz_TwoOuterIonsNeighborTunneling{:.2f}KHz.pdf'.format(self.N, example, 
                                        ions, middle_states ,str(self.equal_distances), self.omega_x/(2*np.pi* 1.e6), 
                                        self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
                        #Save the data:
                        qsave(output1, '{}ions_example{}_DFSions{}_MiddIonsInitState{}_EqualDistanes{}_Radial{:.2f}MHz_Axial{:.2f}KHz_TwoOuterIonsNeighborTunneling{:.2f}KHz'.format(self.N, example, 
                                        ions, middle_states ,str(self.equal_distances), self.omega_x/(2*np.pi* 1.e6), 
                                        self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) 
                                )
                

			 
				  
 			
	
 




    	def DFS_initial_state(self, ions, middle_states=[]):
        	'''Create a pure state where ions[0] and ions[1] are in DFS state and the rest are in 
        	ground state.

        	'''
        	arrays = [ [0] * self.N, [0] * self.N ] 
        	# Note: Don't do arrays = [ [0] * self.N ] * 2, otherwise both arrays[0] and arrays[1] 
        	# will have the same address and changing arrays[0][1] changes arrays[1][0] as well.
 
	        arrays[0][ions[0]-1], arrays[1][ions[1]-1] = 1, 1

	        if middle_states != []: # This works only for ions = [1, N]
	        	for i in range(len(middle_states)):
	        		arrays[0][1+i], arrays[1][1+i]  = middle_states[i], middle_states[i]

                arrays[0][ions[0]-1], arrays[1][ions[1]-1] = 1, 1
        
        	return 1./sqrt(2) * ( ch.ket(arrays[0]) + ch.ket(arrays[1]) )

    	def one_ion_excited_initial_state(self, state, ion):
        	'''Return an initial state with all the ions in ground motional state 
        	except one, given by ion, which is between 1 and self.N.
        	
		'''
		arr = [0] * self.N
  	      	arr[ion-1] = state
        	return ch.ket(arr)
   		
      	def one_ion_entangled_state(self, ion):

    		arr = [[0] * self.N, [0] * self.N]
    		arr[0][ion-1]  = 1
    		arr[1][ion-1] = 0
    		return ( ch.ket(arr[0]) + ch.ket(arr[1]) )/sqrt(2)


    	def self_correlation_op_arr(self, n, ion):
        	''' Return the projection operator of ion number 'ion' to itself 
                for the 'n'-th fock state. 
                ion is an integer between 1 and self.N
                n is an integer smaller than self.M. 
                '''
            	if self.M<n:
                	print "Fock space dimension smaller than phonon number"
		else:
			arr = []
        		for i in range(self.N):
            			if i == ion-1:
               				arr.append(basis(self.M, n) * basis(self.M, n).dag())
            			else:
                			arr.append(qeye(self.M))
        		return arr
		
    	def op_arr(self, ions):
        	arr = []
        	for i in range(self.N):
            		if i == ions[0] - 1:
            			arr.append(basis(self.M, 1)*basis(self.M, 0).dag())
            		elif i == ions[1] - 1:
                		arr.append(basis(self.M, 0)*basis(self.M, 1).dag())
            		else:
                		arr.append(qeye(self.M))
        	return arr

    	def one_ion_1to0_projector_op_arr(self, ion):
        	'''Generate sigmax on the ion. ion is an integer between 1 and self.N'''
        	arr = []
        	for i in range(self.N):
            		if i == ion-1:
                		arr.append(basis(self.M, 0) * basis(self.M, 1).dag())
            		else:
                		arr.append(qeye(self.M))
        	return arr 

    	@classmethod
    	def gen_omegax(cls, u, t12, N):
        	'''Generate omegax under the assumption that the spacings between ions are equal.
        	'''
        	return t12*array(   [ [ u if i == j else 1./(absolute(i-j)**3) for i in range(N)] for j in range(N) ]    )


   	@classmethod
        def generate_couplings(cls, N, omega_x, zposition_arr = [], nearest_neighbor_coupling = 0, 
				axial_freq  = 0, mass = 40 * 1.672621e-27):
        	'''Return the matrix of couplings.
        	Note that nearest_neighbor_coupling is used only when zposition_arr is empty.
        	''' 
        	eps0 = 8.85419e-12
        	echarge = 1.60218e-19
        	if zposition_arr != []:
                	k = echarge**2/(8*mass*omega_x*np.pi*eps0)
                elif nearest_neighbor_coupling != 0:      
                    	k = nearest_neighbor_coupling 
                	zposition_arr = range(N)
                else:
                    	raise Exception("Either ion positions or the nearest neighbor coupling is missing!") 

            	if zposition_arr != []:          
                	t = np.zeros((len(zposition_arr),len(zposition_arr)))
                	for i in range(len(zposition_arr)):
               			for j in range(len(zposition_arr)):
                        		if i != j:
                            			t[i][j] =  k / absolute(zposition_arr[i]-zposition_arr[j]) ** 3
                	return t

        @classmethod
        def generate_local_radial_freqs(cls, omega_x, couplings):
        	"""Return the array of local radial frequencies"""
        	ion_numbers = range(len(couplings[0]))
        	local_radial_freqs = [omega_x for i in ion_numbers]
        	for i in ion_numbers:
          		for j in ion_numbers:
                    		if j != i:
                        		local_radial_freqs[i] -= couplings[i][j]
        	return local_radial_freqs
	
	@classmethod
	def generate_omegax(cls, N, omega_x, zpositions, nearest_neighbor_coupling): 

		couplings = chain.generate_couplings(N, omega_x, zpositions, nearest_neighbor_coupling)
		local_radial_freqs = chain.generate_local_radial_freqs(omega_x, couplings) 
		omegax = np.zeros((N, N))
		#  Construct the matrix of local radial frequencies and couplings 	
		for i in range(N):
			for j in range(N):
				if i == j:
					omegax[i][i] = local_radial_freqs[i]
				else:
					omegax[i][j] = couplings[i][j]
		return omegax
	

	def free_Hamiltonian(self):	
		"""
		Return the free Hamiltonian of phonons in the ion chain
		
		"""    
        	H0 = 0
    		# Write the free Hamiltonian of phonons:
    		for i in range(self.N):
        		for j in range(self.N):
            			H0 += self.omegax[i][j] * (self.a[i]).dag() * self.a[j]
		return H0



	@classmethod
	def eigenvalues(cls, qmat):
		"""Return an array of eigenvalues of the given quantum object"""
		return qmat.eigenstates()[0]	  



	def find_plot_eigs(self, plot):
		"""
		Find eigenvalues, eigenvectors of omegax, save them and plot the results.
		
		""" 
		eigs = LA.eig(self.omegax)
		
		np.savetxt("{}ions_Radial{:.2f}MHz_Axial{:.2f}KHz_ions{}NeighborTunneling{:.2f}KHz_EigenValues.txt".format(self.N,
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), [1,2], self.omegax[0][1]/(2*np.pi* 1.e3) ), eigs[0] )

		np.savetxt("{}ions_Radial{:.2f}MHz_Axial{:.2f}KHz_ions{}NeighborTunneling{:.2f}KHz_EigenVectors.txt".format(self.N,
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), [1,2], self.omegax[0][1]/(2*np.pi* 1.e3) ), eigs[1] )
	

		if plotit:
			eigenfreqs = np.sort([np.real(e) for e in eigs[0]])[::-1] 
	
			plt.figure(11)
		
			plt.xticks(range(1, ch.N+1))
			if equal_distances:
				plt.plot(eigenfreqs/(2*np.pi * 1.e6), 'bs', label='Equal Distances')
			else:
				plt.plot(eigenfreqs/(2*np.pi * 1.e6), 'g^', label='Harmonic Potential')
			plt.xlabel("Normal modes")
			plt.ylabel("Frequencies (MHz)")
			plt.savefig("{}ions_NormalModeFreqs_Radial{:.2f}MHz_Axial{:.2f}KHz_ions{}NeighborTunneling{:.2f}KHz.jpg".format(ch.N,
					 ch.omega_x/(2*np.pi* 1.e6), ch.omega_z/(2*np.pi* 1.e3), [1,2], ch.omegax[0][1]/(2*np.pi* 1.e3) ) )


	def local_thermal_probs(nbar, error = 1.e-2):
	    # Find minimum dimension of effective Hilbert space for a given error and 
	    # return a list with thermal distribution probabilities. 
	    prob_sum = 0
	    nmax = int(nbar)
	    while 1 - absolute(prob_sum) > error:
	        nmax += 1        
	        p_n = [nbar**l/(nbar+1.)**(l+1) for l in range(0,nmax+1)]
	        prob_sum = sum(p_n)
                    
	    return p_n

        def generate_middle_states():
                nmax = 3
                middle_ions_states = [] 
                for nmiddle in range(0, nmax + 1):
                        for i in range(nmax):
                                for j in range(nmax):
                                        if i<=j:
                                                middle_ions_states.append([i, nmiddle, j])
                return middle_ions_states
"""
    def eta(self):
            theta = self.wavevector_angle
            mass =  self.amumass * units.amu
            k = 2.*np.pi/self.laser_wavelength
            eta = k*(units.hbar/(2*mass*2*np.pi*self.trap_frequency))**.5 * np.abs(np.cos(theta*2.*np.pi / 360.0))
            eta = eta.inBaseUnits().value
            return eta



        prob_list = reshape( kron(kron(arr,arr), arr), (1,len(arr)) )
        def global_thermal_probs(error = 5.e-2):
                # Find minimum dimension of effective Hilbert space for a given error and 
                # return a list with thermal distribution probabilities. 
                prob_sum = 0
                nmax = 1
                while (1 - absolute(prob_sum)) > error and nmax<len(prob_list):        
                        p_n = [prob_list[l] for l in range(0,nmax)]
                        prob_sum = sum(p_n)
                        nmax += 1
                return p_n


probs = global_thermal_probs()

[probs[i] for i in range(len(probs)) if probs[i] != probs[i-1]]

#For obtaining the middle ions states for thermal state:
nmax = 3
             ...: middle_ions_states = [] 
     ...: for nmiddle in range(0, nmax + 1):
     ...:     for i in range(nmax):
     ...:         for j in range(nmax):
     ...:             if i<=j:
     ...:                 middle_ions_states.append([i, nmiddle, j])
     ...:                 



"""











if __name__ == '__main__':
	#posit = [] # Position of ion0, ion1, ... on the trap z axis
	nearest_neighbor_coupling =  2 * np.pi * 10.e3
	#omegax = [[2,1],[1,2]] # An example
	# omegax[0][0] = # Initialize local site frequencies x radial direction 
	# omegax[0][1] = # Initialize tunnelings x radial direction
	#u = 2.25e3/6.7
	#omega_x = 2.25e6
	omega_x = 2*np.pi * 2.00e6 #Radial trap frequency with one ion
	omega_z = 2*np.pi * 177.34e3 # 188.10e3 #144.06e3 # 144.04e3
	N = 5 # Number of ions
	eta = .05 # Lamb-Dicke factor
	#M = 3 # Dimension of local Fock space

   	ion = 1 # for example 2,3,4 
	target_ions = [1, N]
	# omegax: Local site frequency (when i = j) and tunnelings (when i != j). Must be real, symmetric.
	
	p = simulation_parameters.simulation_parameters()
	#for middle_ions_fock_state in [[0,0,0]]:
        input_energies = [1,2,3]
    	fock_states_energy_sorted = [ [1, 0,0,0,0], [1, 1,0,0,0], [1, 0,1,0,0], [1, 0,0,1,0], [1, 0,0,0,1], [1, 1,1,0,0], 
                                    [1, 1,0,1,0], [1, 1,0,0,1],[1, 0,1,1,0], [1, 0,1,0,1], [1, 0,0,1,1], 
                                    [1, 2,0,0,0],[1, 0,2,0,0],[1, 0,0,2,0],[1, 0,0,0,2], 
                                    [1, 0,1,1,1], 
                                    [1, 1,0,1,1], [1, 1,1,0,1], [1, 1,1,1,0], 
                                    [1,3,0,0,0], [1,0,3,0,0],[1,0,0,3,0],[1,0,0,0,3],  
                                    [1, 0,0,2,2], [1, 0,2,0,2], [1, 0,2,2,0], [1, 2,0,2,0], [1, 2,2,0,0], [1, 2,0,0,2], 
                                    [1, 1,1,1,1],[1, 2,1,1,1], [1, 1,2,1,1], 
                                    [1, 1,1,2,1],[1, 1,1,1,2], 
                                    
                                    [1, 2,2,1,1], [1, 2,1,2,1], [1, 2,1,1,2], [1, 1, 1,2,2], [1, 1,2,1,2], [1, 2,1,2,2], 
                                    [1, 1,2,2,2], [1, 2,2,1,2], [1, 2,2,2,1], [1, 2,2,2,2], [1, 2,2,0,2], [1, 2,0,2,2], 
                                    [1, 0,2,2,2], [1, 3,2,2,2], [1, 2,3,2,2], [1, 2,2,3,2], [1, 2,2,2,3], [1, 3,3,2,2], 
                                    [1, 3,2,3,2], [1, 3,2,2,3], [1, 2,3,3,2], [1, 2,3,2,3], [1, 2,2,3,3], [1, 3,3,3,2], 
                                    [1, 3,3,2,3], [1, 3,2,3,3], [1, 2,3,3,3], [1, 3,3,3,3], ]  

	for ions_init_fock_state in [0,1,0,0,0]:
        	for equal_distances in [False]:
			if equal_distances:
				zpositions = []
				
			else:
    				zpositions = equilibrium_positions.get_positions(N, omega_z, p)
			
			#  Find couplings and local radial frequencies and put all the values into omegax:
   			omegax = chain.generate_omegax(N, omega_x, zpositions, nearest_neighbor_coupling)     
		
		
			time_scale = 1000.e-6 #100./omegax[0][1]
			times = linspace(0, time_scale, 1000)
			example = 4

			if sum(ions_init_fock_state) <= 2:
				M = 4			
			elif sum(ions_init_fock_state) <= 4:
				M =  3*sum(ions_init_fock_state) # Dimension of local Fock space
			else:
			     	M = 3*sum(ions_init_fock_state) + 5 # Dimension of local Fock space
 	
			#  Set AC-Stark Hamiltonian parameters:  
			# Note: We need around 10KHz rf for carrier
			#rf = 0 #2*np.pi * 800.e3 #2*np.pi * 100.e3  #2*np.pi * 1./500.e-6 # Rabi frequency
			#DELTA = 0 #2*np.pi * 200.e3	# Detuning	
			#omega0 = omegax[ion-1][ion-1]

			# Correction due to AC-Stark on site ion:
			#omegax[ion-1][ion-1] += rf**2/(2*DELTA) * (eta*DELTA)**2/(DELTA**2 - omega0**2) 
		
			# Create a chain:
			ch = chain(omegax, N, M, omega_x, omega_z, equal_distances)	
		    #    middle_ions_fock_state = ions_init_fock_state[1:N-1]
			ch.simulate(target_ions, times, example, ions_init_fock_state, ion)
		
	
	
"""
https://www.youtube.com/watch?v=ZM2lfekLmBE
To be done:
DONE 1. Run all simulations once for equal distances
DONE 2. Import distances from Mike's
DONE 3. Run all simulations once again for harmonic potential mike's distances

DONE Add omega_Z 

Add generate chain methods overloaded with different inputs
Add create trap, with which you set axial_freq and omega_X in the simplest version

Ion trap tool box:
    create trap
    set up trap parameters
    create ion chain
    generate initial state for the chain
    construct operations: such as time evolution, etc
    construct observables   

Reproduce 5 ions with a displacement operator?

Reproduce our actual experiment

DONE Add AC Stark Hamiltonian

Add a graphing function, with capability of labeling axis and setting limits

Add examples as funsions in __main__ not inside the simulation class.
"""   
