"""
An ODE based model of miRNA and Dicer interaction.

Background and data obtained from Tsutsumi et al. 2011, Nat Struct Mol Biol, 10.1038/nsmb.2125

"""

#Required libraries
import numpy as np
import utils
import params
import tqdm
from scipy.integrate import solve_ivp

#import parameters
dt = params.dt
minutes = params.minutes

WT_init = params.WT_init
short_init = params.short_init
dicer_init = params.dicer_init
mirna_init = params.mirna_init
WT_dicer_init = params.WT_dicer_init
short_dicer_init = params.short_dicer_init

init_values = [WT_init, short_init, dicer_init, dicer_init, WT_dicer_init,
                      short_dicer_init, mirna_init, mirna_init]

k1 = params.k1
k_1 = params.k_1
k2 = params.k2
k_2 = params.k_2
k3 = params.k3

theta = np.log([k1, k2, k3])

ka = np.log(k1)
kb = np.log(k2)
kc = np.log(k3)

#functions
def ODE_model(t, init_values, ka, kb, kc):
    """
    This is a function to be passed to an ODE solver.
    
    Args
    t (ndarray):           Time points
    init_values (1darray): Values from function at t
    
    Returns
    WT, short (ndarrays):  Arrays containing concentrations of WT, short miRNA
    """
    WT0, short0, dicer10, dicer20, WT_dicer0, short_dicer0, mirna10, mirna20 = init_values
    k1 = np.exp(ka)
    k2 = np.exp(kb)
    k3 = np.exp(kc)
    
    #WT mirna loop size
    WT = WT_dicer0*k_1 - WT0*dicer10*k1
    WT_dicer = WT0*dicer10*k1 - WT_dicer0*(k_1 + k3)
    dicer1 = WT_dicer0*(k_1 + k3) - WT0*dicer10*k1
    mirna1 = WT_dicer0*k3
    
    #short loop mirna
    short = short_dicer0*k_2 - short0*dicer20*k2
    short_dicer = short0*dicer20*k2 - short_dicer0*(k_2 + k3)
    dicer2 = short_dicer0*(k_2 + k3) - short0*dicer20*k2
    mirna2 = short_dicer0*k3
    
    return WT, short, dicer1, dicer2, WT_dicer, short_dicer,  mirna1, mirna2
    

def conc_change(theta):
    """
    The main model function. Fills out numpy arrays through solving ODEs with the Euler method.
    
    NB! This method bugs out due to numbers too small to handle by computer.
    NB2: Now fixed with v. small stepsize, passing log of theta and transforming when unpacking.
    
    Args
    theta (1darray): Initial reaction rates for each species
    """
    
    k1, k2, k3 = np.exp(theta)
    
    arrays = utils.generate_arrays(species = ['WT', 'dicer1', 'short', 'dicer2', 'WT_dicer', 'short_dicer', 'mirna1', 'mirna2'], 
                                   init_conc = [WT_init, dicer_init, short_init, dicer_init, 0, 0, mirna_init, mirna_init],
                                   dt = dt, minutes = minutes)
    
    WT = arrays['WT']
    dicer1 = arrays['dicer1']
    short = arrays['short']
    dicer2 = arrays['dicer2']
    WT_dicer = arrays['WT_dicer']
    short_dicer = arrays['short_dicer']
    mirna1 = arrays['mirna1']
    mirna2 = arrays['mirna2']
    
    for i in tqdm.tqdm(range(1, len(WT))):
         #force 0 if concentration negative
            ##does not work
        #WT[i-1] = max(WT[i-1], 0)
        #dicer1[i-1] = max(dicer1[i-1], 0)
        #WT_dicer[i-1] = max(WT_dicer[i-1], 0)
        #mirna1[i-1] = max(mirna1[i-1], 0)
        
        #short[i-1] = max(short[i-1], 0)
        #dicer2[i-1] = max(dicer2[i-1], 0)
        #short_dicer[i-1] = max(short_dicer[i-1], 0)
        #mirna2[i-1] = max(mirna2[i-1], 0)
        
        #wild type mirna loop
        WT[i] = WT[i-1] + dt*(WT_dicer[i-1]*k_1 - WT[i-1]*dicer1[i-1]*k1)
        dicer1[i] = dicer1[i-1] + dt*(WT_dicer[i-1]*(k_1 + k3) - WT[i-1]*dicer1[i-1]*k1)
        WT_dicer[i] = WT_dicer[i-1] + dt*(WT[i-1]*dicer1[i-1]*k1 - WT_dicer[i-1]*(k_1 + k3))
        mirna1[i] = mirna1[i-1] + dt*(WT_dicer[i-1]*k3)
        
        #short mirna loop
        short[i] = short[i-1] + dt*(short_dicer[i-1]*k_2 - short[i-1]*dicer2[i-1]*k2)
        dicer2[i] = dicer2[i-1] + dt*(short_dicer[i-1]*(k_2 + k3) - short[i-1]*dicer2[i-1]*k2)
        short_dicer[i] = short_dicer[i-1] + dt*(short[i-1]*dicer2[i-1]*k2 - short_dicer[i-1]*(k_2 + k3))
        mirna2[i] = mirna2[i-1] +dt*(short_dicer[i-1]*k3)
        
    return WT, short


def frac_diced(theta):
    """
    Function to calculate the fraction of WT, short miRNA that has been diced from stepwise function
    """
    
    WT, short = conc_change(theta)
    
    arrays = utils.generate_arrays(species = ['WT_diced', 'short_diced'], init_conc = [0, 0],
                                  dt = dt, minutes = minutes)
    WT_diced = arrays['WT_diced']
    short_diced = arrays['short_diced']
    
    for i in range(len(WT)):
        WT_diced[i] = (WT[0] - WT[i]) / WT[0]
        short_diced[i] = (short[0] - short[i]) / short[0]
        
    return WT_diced, short_diced

def frac_diced_ODE(theta):
    """
    Function to calculate the fraction of WT, short miRNA that has been diced from ODE function
    """
    k1, k2, k3 = theta
    
    sol = solve_ivp(ODE_model, (0, int(params.minutes)), init_values, args = (k1, k2, k3))
    
    WT, short, dicer1, dicer2, WT_dicer, short_dicer,  mirna1, mirna2 = sol.y
    
    arrays = utils.generate_arrays(species = ['WT_diced', 'short_diced'], init_conc = [0, 0],
                                  length = [int(len(WT)), int(len(short))])
    WT_diced = arrays['WT_diced']
    short_diced = arrays['short_diced']
    
    for i in range(len(WT)):
        WT_diced[i] = (WT[0] - WT[i]) / WT[0]
        short_diced[i] = (short[0] - short[i]) / short[0]
        
    return WT_diced, short_diced, sol.t