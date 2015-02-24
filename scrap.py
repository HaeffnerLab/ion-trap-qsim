#!/usr/env/ python





def get_H0(omegax_interactionPic, rel_detuning):   
    #rel_detuning is the distance of laser frequency relative to the sideband, positive is less than sideband
    s   =   0
    for i in range(N):
        for j in range(N):
            s  +=   omegax_interactionPic[i][j] * a[i].dag()*a[j]
    return rel_detuning *  sum([a[j].dag()*a[j] for j in range(N)]) + s
        
    
def get_H(ion_num, phase):
    return  1.j * eta*omega * ( exp(-1.j*phase) * sigma_plus[ion_num] * a[ion_num].dag() - exp(1.j*phase) * sigma_minus[ion_num] * a[ion_num] ) 



def do_lattice_with_gen():
    exp_final  =  np.zeros(len(times_free))
    ion_num    =  0
    rel_detuning =  DELTA
    H0        =  get_free_hamiltonian() 
    #pulse pair 1#
    U1        =  propagator_Trotter(H0 + get_H(0, 0), times_pulse, 'final')[-1] 
    #pulse pair 2#
    U2        =  propagator_Trotter(H0 + get_H(0, second_pair_pulses_phase), times_pulse, 'final')[-1]  
    rho0      =  U1 * rho0 * U1.dag()
    
    gen       =  propagator_gen( H0, time_precision )
    for j,t in enumerate(times_free):
        #U0       =  propagator_Trotter( H0, times_free )
        U0_t      =  gen.next()
        u         =  U2 * U0_t 
        
        exp_final[j]    = real(expect ( excited_state_pop[ion_num], u * rho0 * u.dag() ) ) 
        
    #rho_final = rho_final/laser_detuning_mc_num
    #exp_final  = exp_final/laser_detuning_mc_num
    #return rho_final
    return exp_final



def do_lattice_with_gen(self, times_free, ):
    exp_final  =  np.zeros(len(times_free))
    ion_num    =  0
    rel_detuning =  DELTA
    H0        =  get_H0(omegax, rel_detuning) 
    #pulse pair 1#
    U1        =  propagator_Trotter(H0 + get_H(0, 0), times_pulse, 'final')[-1] 
    #pulse pair 2#
    U2        =  propagator_Trotter(H0 + get_H(0, second_pair_pulses_phase), times_pulse, 'final')[-1]  
    rho0      =  U1 * rho0 * U1.dag()
    
    gen       =  propagator_gen( H0, time_precision )
    for j,t in enumerate(times_free):
        #U0       =  propagator_Trotter( H0, times_free )
        U0_t      =  gen.next()
        u         =  U2 * U0_t 
        
        exp_final[j]    = real(expect ( excited_state_pop[ion_num], u * rho0 * u.dag() ) ) 
        
    #rho_final = rho_final/laser_detuning_mc_num
    #exp_final  = exp_final/laser_detuning_mc_num
    #return rho_final
    return exp_final

ion_num  =  0
DELTA    =  omegax_init - sum( [omegax[i][i] for i in range(N)] )/N #Laser freq positive from the highest sideband line, i.e., radial frequency.  

#exp_final = do_1param_mc_with_gen()
#exp_final = do_complete_2param_mc()
lattice_exp_final  = do_lattice_with_gen() 






