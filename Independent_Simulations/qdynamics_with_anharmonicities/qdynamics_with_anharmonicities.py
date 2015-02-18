#!/usr/bin/env python

import matplotlib.pylab as plt
from qutip import *
from scipy import *
import numpy as np
from numpy import *
from equilibrium_positions import equilibrium_positions
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
		# Create local Fock space destruction operators:
	        self.a = [self.ann_op(i) for i in range(self.N)]
	
		
	def simulate(self, ions, times, example, ion = 1, OMEGA=0, DELTA=0, delta_radial_freq=0, 
                        omegax_gradient_at_100microns=2*np.pi * 10.e3):
	        N = self.N
        	M = self.M
            	self.delta_radial_freq = delta_radial_freq
        	self.omegax_gradient_at_100microns  =  omegax_gradient_at_100microns
                
		if example == 1: #initial state (|01>+|10>)/sqrt(2), op = |01><10| + |10><01|, for ions[0] and ions[1]   
			psi0 = self.DFS_initial_state(ions)
        		op_term1 = tensor(self.op_arr(ions))
			op = op_term1 + op_term1.dag()
			# Find time evolution of psi0 and op with free Hamiltonian
    			H = self.free_Hamiltonian()
    			output1 = mcsolve(H, psi0, times, [], [op])
		

    		elif example == 2: 
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
			# Find time evolution of psi0 and op with free Hamiltonian
    			H = self.free_Hamiltonian()
    			output1 = mcsolve(H, psi0, times, [], [op])
		

        	elif example == 4:
			psi0 = self.one_ion_excited_initial_state(ion)
    			op_term1 = tensor(self.self_correlation_op_arr(ion))
    			op = op_term1
			AC_Stark_ion = ion
			H = self.free_Hamiltonian() #+ self.AC_Stark_Hamiltonian(AC_Stark_ion, DELTA, OMEGA)			
	                
			# Find time evolution of psi0 and op
    			output1 = mcsolve(H, psi0, times, [], [op])
			
		elif example == 5:
                        #Correct only the radial freq of the two end ions by self.delta_radial_freq

                        psi0 = self.one_ion_excited_initial_state(ion)
                        op_term1 = tensor(self.self_correlation_op_arr(ion))
                        op = op_term1
                        #modify local radial frequency of last two ions by self.delta_radial_freq:
                        self.omegax[0][0] +=  self.delta_radial_freq
                        self.omegax[self.N-1][self.N-1] += self.delta_radial_freq

                        H = self.free_Hamiltonian() #+ self.AC_Stark_Hamiltonian(AC_Stark_ion, DELTA, OMEGA)                    
                        
                        # Find time evolution of psi0 and op
                        output1 = mcsolve(H, psi0, times, [], [op])		


            	elif example == 6:
                        #Correct only the radial freq of all ions

                        psi0 = self.one_ion_excited_initial_state(ion)
                        op_term1 = tensor(self.self_correlation_op_arr(ion))
                        op = op_term1
                        #modify local radial frequency of all ions:
                        omegax_corrections = self.get_radial_freq_corrections()
                        for i in range(self.N):
                            self.omegax[i][i] +=  omegax_corrections[i]

                        H = self.free_Hamiltonian()               
                        
                        # Find time evolution of psi0 and op
                        output1 = mcsolve(H, psi0, times, [], [op])     

    			

    		#Plot and save simulation data:    
    		if example == 1: 
			plt.figure(example*self.N)
			plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
	        	plt.xlabel('Time ($\mu$s)')
			plt.ylabel('Expectation value of op')
	        	plt.ylim(-1,1)
			plt.title('Expectation value of op = |01><10|+|10><01| for ions {} in a chain of {} ions, with max = {}'.format(ions,					N, max(output1.expect[0])))
			plt.xlim(0, times[-1]*1.e6)
			#Save the graph:
			plt.savefig('{}ions_ions{}_example{}_EqualDistanes{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg'.format(self.N, ions, 
                example,  str(self.equal_distances), self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[ions[0]-1][ions[1]-1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
			#Save the data:
			qsave(output1, '{}ions_ions{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz_data'.format(self.N, ions, example,  str(self.equal_distances), self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[ions[0]-1][ions[1]-1]/(2*np.pi* 1.e3) ) )
		
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
			plt.ylabel('Probability')
 	        	plt.ylim(0,1)
			plt.title('Probability of finding ion number {} in excited state, in a chain of {} ions, with max = {}'.format(ion, 
				N, max(output1.expect[0])))
        		plt.xlim(0, times[-1]*1.e6)
			#Save the graph:
			plt.savefig('data/{}ions_ion{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg'.format(self.N, ion, example,  str(self.equal_distances),
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
			qsave(output1, 'data/{}ions_ion{}_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz_data'.format(self.N, ion, example,  str(self.equal_distances),
					 self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )
 
   	        elif example == 5:
                        plt.figure(example*(self.N+25)) # 25 is the max number of runs
                        plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
                        plt.xlabel('Time ($\mu$s)')
                        plt.ylabel('Probability')
                        plt.ylim(0,1)
                        plt.title('Probability of finding ion number {} in excited state, in a chain of {} ions, with max = {}'.format(ion, 
                                N, max(output1.expect[0])))
                        plt.xlim(0, times[-1]*1.e6)
                        #Save the graph:
                        plt.savefig('data/{}ions_radial_freq_modification{:.2f}KHz_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg'.format(self.N, 1.e-6*self.delta_radial_freq/(2*np.pi* 1.e3),  example,  str(self.equal_distances),
                                         self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
                        qsave(output1, 'data/{}ions_radial_freq_modification{:.2f}KHz_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz_data'.format(self.N, 1.e-6*self.delta_radial_freq/(2*np.pi* 1.e3), example,  str(self.equal_distances),
                                         self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )

            	elif example == 6:
                        plt.figure(example*(self.N+25)) # 25 is the max number of runs
                        plt.plot([x*1.e6 for x in output1.times], output1.expect[0])
                        plt.xlabel('Time ($\mu$s)')
                        plt.ylabel('Probability')
                        plt.ylim(0,1)
                        plt.title('Probability of finding ion number {} in excited state, in a chain of {} ions, with max = {}'.format(ion, 
                                N, max(output1.expect[0])))
                        plt.xlim(0, times[-1]*1.e6)
                        #Save the graph:
                        plt.savefig('data/{}ions_all_ions_radial_freq_modification_with_gradient_at_100microns_of_{:.2f}KHzPerMicrons_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg'.format(self.N, 1.e-6*self.omegax_gradient_at_100microns/(2*np.pi* 1.e3), example,  str(self.equal_distances),
                                         self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ), bbox_inches='tight')
                        qsave(output1, 'data/{}ions_all_ions_radial_freq_modification_with_gradient_at_100microns_of_{:.2f}KHzPerMicrons_example{}_EqualDistances{}_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz_data'.format(self.N, 1.e-6*self.omegax_gradient_at_100microns/(2*np.pi* 1.e3), example,  str(self.equal_distances),
                                         self.omega_x/(2*np.pi* 1.e6), self.omega_z/(2*np.pi* 1.e3), self.omegax[0][1]/(2*np.pi* 1.e3) ) )

                

       
		#plt.show()

    	def DFS_initial_state(self, ions):
        	'''Create a pure state where ions[0] and ions[1] are in DFS state and the rest are in 
        	ground state.

        	'''
        	arrays = [ [0] * self.N, [0] * self.N ] 
        	# Note: Don't do arrays = [ [0] * self.N ] * 2, otherwise both arrays[0] and arrays[1] 
        	# will have the same address and changing arrays[0][1] changes arrays[1][0] as well.
 
	        arrays[0][ions[0]-1], arrays[1][ions[1]-1] = 1, 1
        	return 1./sqrt(2) * ( ch.ket(arrays[0]) + ch.ket(arrays[1]) )

    	def one_ion_excited_initial_state(self, ion):
        	'''Return an initial state with all the ions in ground state 
        	except one, given by ion, which is between 1 and self.N.
        	'''
        	arr = [0] * self.N
        	arr[ion-1] = 1
        	return ch.ket(arr)
   
      	def one_ion_entangled_state(self, ion):

    		arr = [[0] * self.N, [0] * self.N]
    		arr[0][ion-1]  = 1
    		arr[1][ion-1] = 0
    		return ( ch.ket(arr[0]) + ch.ket(arr[1]) )/sqrt(2)


    	def self_correlation_op_arr(self, ion):
        	'''Generate sigmax on the ion. ion is an integer between 1 and self.N'''
        	arr = []
        	for i in range(self.N):
            		if i == ion-1:
               			arr.append(basis(self.M, 1) * basis(self.M, 1).dag())
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

	"""
	def evolve( self, H, t, psi0 = self.ket([0]*self.N) ):
		
		Evolve the initial state for the given H and duration (t) 
		and return the final state.
		
		times = linspace(0, t, 1000) 
		output = mcsolve(H, psi0, times, [], [])
		
		return output.states[-1]
	"""
	
	def AC_Stark_Hamiltonian( self, AC_Stark_ion, DELTA, OMEGA ):
		"""Evolve the initial state for the given duration (t), detuning(DELTA), rabi frequency (OMEGA)
		and return the final state."""

		arr1 = [0]*self.N
		arr1[AC_Stark_ion - 1] = 1
		arr0 = [0]*self.N

		ket1 = self.ket(arr1)
		ket0 = self.ket(arr0)
		H_ACStark = -OMEGA**2/(4*DELTA) * ( ket1*ket1.dag() - ket0*ket0.dag()  ) 	

		return H_ACStark

	
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


    	def get_radial_freq_corrections(self):
        	p = simulation_parameters.simulation_parameters()
        	zpositions = equilibrium_positions.get_positions(self.N, self.omega_z, p)
        	self.delta_omegax = self.omegax_gradient_at_100microns * zpositions**2./(2*100.e-6)

        	return self.delta_omegax 



if __name__ == '__main__':
	#posit = [] # Position of ion0, ion1, ... on the trap z axis
	nearest_neighbor_coupling =  2 * np.pi * 10.e3
	#omegax = [[2,1],[1,2]] # An example
	# omegax[0][0] = # Initialize local site frequencies x radial direction 
	# omegax[0][1] = # Initialize tunnelings x radial direction
	#u = 2.25e3/6.7
	#omega_x = 2.25e6
	omega_x = 2*np.pi * 2.00e6 #Radial trap frequency with one ion
	omega_z = 2*np.pi * 100.00e3
        #omegax_gradient_at_100microns_arr = arange(0.2,0, -.05) * 2*np.pi * 10.e3/(1.e-6)  # Per microns 
    	omegax_gradient_at_100microns_arr = np.array([1., 0.2, .1, 0.05, 0.03, 0])
	N = 10 # Number of ions
	
    	M = 2 # Dimension of local Fock space
   	ion = 1
	target_ions = [1, 10]
	# omegax: Local site frequency (when i = j) and tunnelings (when i != j). Must be real, symmetric:
	# omegax = t12*array([[u, 1],[1, v]]) # 2 ions
	#omegax = t12*array([[u, 1, 1./8, 1./27],[1, v, 1., 1./8],[1./8, 1, u, 1.], [1./27, 1./8, 1., u]]) # 4 ions or less
	
	p = simulation_parameters.simulation_parameters()
    
    	for omegax_gradient_at_100microns in omegax_gradient_at_100microns_arr:	
    		for equal_distances in [False]:
    		
    			if equal_distances:
    				zpositions = []
    			
    			else:
        			zpositions = equilibrium_positions.get_positions(N, omega_z, p)
    			
    			#  Find couplings and local radial frequencies and put all the values into omegax:
       			omegax = chain.generate_omegax(N, omega_x, zpositions, nearest_neighbor_coupling)     
    		
    			ch = chain(omegax, N, M, omega_x, omega_z, equal_distances)	
    			time_scale = 0.7*100./omegax[0][1]
    			times = linspace(0, time_scale, 1000)
    			example = 6
    				
    			#  Set AC-Stark Hamiltonian parameters:  
    			# Note: We need around 10KHz rf for carrier
    			rf = 2*np.pi * 100.e3  #2*np.pi * 1./500.e-6
    			DELTA = 0 #20.* omegax[0][1]		
    			delta_radial_freq = +2*np.pi * 100.e3

    			ch.simulate(target_ions, times, example, ion, rf, DELTA, delta_radial_freq, omegax_gradient_at_100microns)

    			#eigs = LA.eig(ch.omegax)
    			#f = open("3ions_Axial144KHz_Radial2.25MHz_[1,2]ionsTunneling10KHz.txt", 'w')
    			#f.write(str(eigs))
    			#f.close()
    			#np.savetxt("{}ions_Radial{:.2f}MHz_Axial{:.2f}KHz_ions{}NeighborTunneling{:.2f}KHz_EigenValues.txt".format(ch.N,
    			#			 ch.omega_x/(2*np.pi* 1.e6), ch.omega_z/(2*np.pi* 1.e3), target_ions, ch.omegax[0][1]/(2*np.pi* 1.e3) ), eigs[0] )

    			#np.savetxt("{}ions_Radial{:.2f}MHz_Axial{:.2f}KHz_ions{}NeighborTunneling{:.2f}KHz_EigenVectors.txt".format(ch.N,
    			#			 ch.omega_x/(2*np.pi* 1.e6), ch.omega_z/(2*np.pi* 1.e3), target_ions, ch.omegax[0][1]/(2*np.pi* 1.e3) ), eigs[1] )


    			"""
    			eigenfreqs = np.sort([np.real(e) for e in eigs[0]])[::-1] 
    	
    			plt.figure(11)
    		
    			plt.xticks(range(1, ch.N+1))
    			if equal_distances:
    			plt.plot(eigenfreqs/(2*np.pi * 1.e6), 'bs', label='Equal Distances')
    			else:
    			plt.plot(eigenfreqs/(2*np.pi * 1.e6), 'g^', label='Harmonic Potential')
    			plt.xlabel("Normal modes")
    			plt.ylabel("Frequencies (MHz)")
    			plt.savefig("{}ions_NormalModeFreqs_Radial{:.2f}MHz_Axial{:.2f}KHz_NeighborTunneling{:.2f}KHz.jpg".format(ch.N,
    					 ch.omega_x/(2*np.pi* 1.e6), ch.omega_z/(2*np.pi* 1.e3), ch.omegax[0][1]/(2*np.pi* 1.e3) ) )
    			"""
			#	plt.legend([p1,p2], ["Equal distances", "Harmonic potential"])
	
	
'''https://www.youtube.com/watch?v=ZM2lfekLmBE
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
'''   
