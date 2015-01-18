#!/usr/bin/env python
import parameters as p
from ion_trap import IonTrap
from ions import Ion, Chain
from experiment import Experiment
from observables import Observable


class ARamseyExperiment(Experiment):

    def __init__(self, trap, chain, lasers, pulse_sequence, observables,
                 data_storage_path = '/home/trxw/Documents/dfr/codes/quantum_play/ChainSimulator/Data_Storage' ):
       
        self.trap               = trap
        self.chain              = chain
        self.lasers             = lasers
        self.pulse_sequence     = pulse_sequence
        self.observables        = observables                        
        self.data_storage_path  = data_storage_path 

        #super(ARamseyExperiment, self).__init__()

    


if __name__ == '__main__':

        num_of_ions             =  3
        motional_hilbert_space_dim  =  8
        chain                   =  Chain( num_of_ions, motional_hilbert_space_dim )

        #Define trap:
        dummy_trap              =  IonTrap( p.trap_parameters )
        
        dummy_trap.load( chain )

        #Setup lasers:
        lasers                  =  {laser1                  =  Laser(ion_number = 1, freq = -1, intensity = .1)
                                    laser2                  =  Laser(ion_number = 5, freq = 1, intensity = .1)
                                    }

        pulse_sequence          =  { (laser1, laser1_time), 
                                     (laser2, laser2_time),
                                     (laser1, laser1_time),
                                     (laser2, laser2_time),
                                    }

        observables             =  { Observable('ion1', 'electronic_state')

                                    }



        experiment_parameters   =  {'trap'               : trap, 
                                    'chain'              : chain,
                                    'lasers'             : lasers,
                                    'pulse_sequence'     : pulse_sequence,
                                    'observables'        : observables,                         
                                    'data_storage_path'  : '/home/trxw/Documents/dfr/codes/quantum_play/ChainSimulator/Data_Storage'
                                    }

        
        experiment              =  ARamseyExperiment( **experiment_parameters )



        #Initialize the chian quantum states:


        #Run the experiment and collect data:
        experiment.run()

        #Plot
        experiment.plot(observables[0])

