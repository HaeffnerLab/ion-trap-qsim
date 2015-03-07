from __future__ import division, absolute_import, print_function, unicode_literals
from qutip import *

class OperatorZoo(object):

    def __init__(self):


        M             =   self.chain.ion_motional_hilbert_space_dim
        

        self.a_func        = lambda i: tensor( qeye(2), tensor( [ destroy(M) if j == i else qeye(M) for j in range(self.chain.num_of_ions) ] ) )
                                     
        

        self.a                =  [self.a_func(j) for j in range(self.chain.num_of_ions)]

        #Normal modes expressed in terms of local mode a's:
        self.D_func        = lambda i: sum( [ self.chain.normal_in_local_modes[i][j] * self.a[j] for j in range(self.chain.num_of_ions) ] )
         
        self.D             = [ self.D_func(i) for i in range(self.chain.num_of_ions) ]

        


        self.excited_state_pop_func  =  lambda i: tensor( create(2) * destroy(2),
                                                 tensor([ qeye(M) for j in range(self.chain.num_of_ions)] )                                       
                                               )
        self.excited_state_pop       =  [self.excited_state_pop_func(j) for j in range(self.chain.num_of_ions)]



        self.sigma_plus_func  =  lambda i:  tensor( create(2),
                                           tensor( [ qeye(M) for j in range(self.chain.num_of_ions)] ) 
                                        )
        self.sigma_plus       =  [ self.sigma_plus_func(j) for j in range(self.chain.num_of_ions) ]

        self.sigma_minus_func  =  lambda i:  tensor( destroy(2),
                                           tensor( [ qeye(M) for j in range(self.chain.num_of_ions)] ) 
                                        )
        self.sigma_minus      =  [self.sigma_minus_func(j) for j in range(self.chain.num_of_ions)]

        #n        =  lambda i: a(i).dag() * a(i) # local phonon number, N>i>=0

        self.motion_identity_op =  tensor( [qeye(M) for j in range(self.chain.num_of_ions)] )
                               

        #self.sigmaz          =    tensor( sigmaz ,
        #                                self.motion_identity_op    

        #                                )
        self.sigmaz_func  =  lambda i:  tensor( sigmaz,
                                           tensor( [ qeye(M) for j in range(self.chain.num_of_ions)] ) 
                                        )
        self.sigmaz      =  [self.sigma_minus_func(j) for j in range(self.chain.num_of_ions)]



    @property
    def spin_identity_operators_arr(self):
        pass



    def ket(self, arr): 
        return tensor([basis(self.M, e) for e in arr])

    def get_destruction(self, element_num):
        """ Generates a destruction operator which acts on the self.hilbert_space_dim_array's element_num
        element subspace of system Hilbert space. 
        element_num  starts from 1. 

        """
        if element_num > self.chain.num_of_ions:
            print("In get_destruction element_num must be between 1 and number of ions in the chain")
            return None
        else:     
            arr  =   []
            for i, dim in enumerate(self.hilbert_space_dim_array):
                arr   +=   [ qeye(dim) ]


    def get_creation(self, element_num):

        """ Generates a destruction operator which acts on the self.hilbert_space_dim_array's element_num
        element subspace of system Hilbert space. 
        element_num  starts from 1. 

        """
        if element_num > self.chain.num_of_ions:
            print("In get_destruction element_num must be between 1 and number of ions in the chain POTATO")
            return None
        else:     
            arr  =   []
            for i, dim in enumerate(self.hilbert_space_dim_array):
                arr   +=   [ qeye(dim) ]

    @property
    def chain_electronic_states_identity_op(self):
        pass