
import pickle
from numpy import *
from scipy.optimize import basinhopping


import os,sys,inspect


from scipy import linalg as LA 
from scipy import optimize
import numpy as np
from simulation_parameters import SimulationParameters, OptimizationParameters
from numpy.random import random
from ions import Chain
from ion_trap import *
import matplotlib.pylab as plt
from spectral_optimization import *
from scipy.optimize import *
import sys


if __name__ == '__main__':

    N = 15
    chain = Chain(N)
    trap = IonTrap()
    trap.set_trap_radial_freq(2.e6, 'Y')
    trap.set_trap_radial_freq(2.e6, 'X')

    used_electrodes    = range(1,24)
    seed_electrode_voltages = array( [6.938899999999999955e-17, 8.791899999999999715e-02, 4.813600000000000101e-01, 8.047299999999999454e-01,
                                -1.882300000000000084e+00, 8.013900000000000468e-01, 5.024399999999999977e-01, 8.146799999999999875e-02,
                                2.183300000000000171e-02, 0.000000000000000000e+00, 0.000000000000000000e+00, 0.000000000000000000e+00,
                                3.020199999999999968e-02,7.582099999999999951e-02,-2.652300000000000213e-01,-3.169699999999999740e-01,
                                -2.666000000000000036e-01,7.831100000000000561e-02,3.159700000000000009e-02,8.757400000000000254e-03,
                                0.000000000000000000e+00,0.000000000000000000e+00,6.395599999999999896e-02] )

    trap.load(chain)

    op = Optimization(trap)
    func = op.cost_function_voltages_and_max_voltages_exp_punish



    x0 = seed_electrode_voltages

    minimizer_kwargs = {"method" : "BFGS"}


    orignal_stdout = sys.stdout
    f   = file('Optimization/Results/global_optimization_results_2modes_isolation.txt', 'w')
    sys.stdout = f

    ret = basinhopping(func, x0, minimizer_kwargs=minimizer_kwargs, niter=10000)


    f.close()




