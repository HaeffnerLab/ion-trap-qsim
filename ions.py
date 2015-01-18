#!/usr/bin/env python
from qutip import *
from scipy import *
import numpy as np
from error_handling import *

class Ion:
        '''
        Create an Ion object with ion_number (>=0), 
        motional_state (>=0) and electronic_state (either 0 or 1).

        '''
        def __init__( self, ion_number, ion_motional_hilbert_space_dim, ion_electronic_hilbert_space_dim = 2, state_type='pure'):

                self.ion_number                        =  ion_number
                self.ion_motional_hilbert_space_dim    =  ion_motional_hilbert_space_dim
                self.ion_electronic_hilbert_space_dim  =  ion_electronic_hilbert_space_dim
                self.state_type                        =  state_type

        def __initialize_ion_state(func):

                def wrapper(self, state, hilbert_space_dim):

                        if self.state_type == 'pure':
                                if hilbert_space_dim <= state:
                                        raise Dimensionerror( "Ion number %s\n" % self.ion_number )
                                else:
                                        func(self, state, hilbert_space_dim)

                        elif self.state_type == 'density_operator':
                                print "To be coded!"

                return wrapper


        @__initialize_ion_state
        def initialize_ion_motional_state(self, ion_motional_state, motional_hilbert_space_dim):

                self.ion_motional_state     = basis( motional_hilbert_space_dim, ion_motional_state )


        @__initialize_ion_state
        def initialize_ion_electronic_state(self, ion_electronic_state, electronic_hilbert_space_dim):

                self.ion_electronic_state   = basis( electronic_hilbert_space_dim, ion_electronic_state )


        def get_position(self):

                return self.position

        def set_position(self, position):
                self.position   =  position




class Chain:

        def __init__(self, N, ion_motional_hilbert_space_dim, ion_electronic_hilbert_space_dim = 2, state_type = 'pure'):

                self.num_of_ions                       =  N 
                self.ion_motional_hilbert_space_dim    =  ion_motional_hilbert_space_dim
                self.ion_electronic_hilbert_space_dim  =  ion_electronic_hilbert_space_dim
                self.Ions                              =  [Ion(i, self.ion_motional_hilbert_space_dim, self.ion_electronic_hilbert_space_dim, state_type) for i in range(N)]
                

        def set_chain_electronic_states( self, **kwargs): 
                ''' Initialize the initial state of ions given in input,
                Example for input format:  kwargs = {'ion1': 0, 'ion5': 1}
                for a chain of 5 ions. If ion is not specified, it is assumed 
                that it is not involved in the dynamics (So defined lasers 
                should not interact with those states.). Also, if it is not involved 
                in the interactions it must not be included to avoid exceeding memory.
                Note that naming ions start from 1 in the input arguments.
                '''
                for k in kwargs:
                        try:
                                if int(k[k.find('ion')+3:]) - 1 >= 0:
                                        self.Ions[ int(k[k.find('ion')+3:]) - 1  ].initialize_ion_electronic_state( kwargs[k], self.ion_electronic_hilbert_space_dim )
                                else:
                                        "Ion numbering starts at 1."
                        except ValueError, e:
                                print "Error: Ions name formatting must be as 'ion10' ", e


        def set_chain_motional_states( self, *args): 
                ''' Initialize the initial state of ions given in input,
                Example for input format:  args = (0, 0, 0, 1, 0)
                for a chain of 5 ions.If ion is not specified, it is assumed 
                that it is not involved in the dynamics (So defined lasers 
                should not interact with those states.). Note that naming ions start 
                from 1 in the input arguments.
                '''
                if len(args) != self.num_of_ions:
                        print "Number of arguments must be equal to the number of ions in the chain."
                else:
                        for i in range(len(args)):
                                self.Ions[i].initialize_ion_motional_state( args[i], self.ion_motional_hilbert_space_dim )


        def set_zpositions( self, zpositions ):
                if len(zpositions) == self.num_of_ions:
                        for i in range(self.num_of_ions):
                                self.Ions[i].set_position( zpositions[i] )
                else:
                        print "Number of given ions don't match with number of elements in zpositions array."

        def get_positions(self):
                
                return [ion.get_position() for ion in self.Ions]
                