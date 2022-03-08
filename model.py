"""
An ODE based model of miRNA and Dicer interaction.

Background and data obtained from Tsutsumi et al. 2011, Nat Struct Mol Biol, 10.1038/nsmb.2125

"""

#Required libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import cma
import os
import utils

## Model setup
#timesteps
dt = 0.0001
minutes = 60

Kd_wt = 25.4 #nM, experimental values
Kd_short = 147.7 #nM, experimental values

#K_d = k_off / k_on

WT_init = 1 #nM, from experimental setup
short_init = 1 #nM
dicer_init = 5 #nM
mirna = 0

k1 = 5
k_1 = Kd_wt * k1
k2 = 5
k_2 = Kd_short * k2
k3 = 5

theta = [k1, k2, k3]

#data fig 1
time = np.array([0, 5, 10, 20, 40, 60])
WT_y = np.array([0, 0.11144276160503169, 0.16566679779700877, 0.23905143587726366, 0.2954956726986665, 0.2946793863099961])
short_y = np.array([0, 0.0033684107002276975, 0.007599822974028003, 0.010019177812737812, 0.009603658536577298, 0.01242378048779691])

#functions
def conc_change(theta):
    """
    The main model function. Fills out numpy arrays through solving ODEs with the Euler method.
    
    Args
    theta (array of floats): Initial reaction rates for each species
    """
    
    k1, k2, k3 = np.log(theta)
    
    arrays = utils.generate_arrays(species = ['WT', 'dicer1', 'short', 'dicer2', 'WT_dicer', 'short_dicer', 'mirna1', 'mirna2'], 
                                   init_conc = [WT_init, dicer_init, short_init, dicer_init, 0, 0, mirna_init, mirna_init],
                                   dt = dt, minutes = minutes)
    
    WT = arrays['WT']
    dicer1 = arrays['dicer1']
    short = arrays['short']
    dicer2 = arrays['dicer2']
    WT_dicer = arrays9['WT_dicer']
    short_dicer = arrays['short_dicer']
    mirna1 = arrays['mirna1']
    mirna2 = arrays['mirna2']
    
    for i in range(1, len(WT)):
        #wild type mirna loop
        WT[i] = WT[i-1] + dt*( )
        dicer1[i] = dicer[i-1] + dt*( )
        WT_dicer[i] = WT_dicer[i-1] + dt*
        mirna1[i] = mirna1[i-1] + dt*
        
        #short mirna loop
        short[i] = short[i-1] + dt*( )
        dicer2[i] = dicer[i-1] + dt*( )
        short_dicer[i] = shprt_dicer[i-1] + dt*( )
        mirna2[i] = mirna2[i-1] +dt*( )
        
    return WT, short


def frac_diced(theta):
    """
    Function to calculate the fraction of WT, short miRNA that has been diced
    """
    
    WT, short = conc_change(theta)
    
    arrays = utils.generate_arrays(species = ['WT_diced', 'short_diced'], init_conc = [0, 0],
                                  dt = dt, minutes = miinutes)
    WT_diced = arrays['WT_diced']
    short_diced = arrays('short_diced')
    
    for i in range(1, len(WT)):
        