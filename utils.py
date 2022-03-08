"""
A selection of functions required for running and optimisation of the model

"""

import numpy as np
import math

#data fig 1
time = np.array([0, 5, 10, 20, 40, 60])
WT_y = np.array([0, 0.11144276160503169, 0.16566679779700877, 0.23905143587726366, 0.2954956726986665, 0.2946793863099961])
short_y = np.array([0, 0.0033684107002276975, 0.007599822974028003, 0.010019177812737812, 0.009603658536577298, 0.01242378048779691])

data = np.array([[WT_y],
                 [short_y]])

def generate_arrays(*, species, init_conc, dt = 0.01, minutes = 60):
    """
    A function that takes initially provided parameters and produces a
    selection of numpy arrays to be filled through solving model ODEs
    using the Euler method.
    
    Args
    species (array, string):     List of species names to be used to 
    init_conc (array of floats): List of inital concentrations in nM
    dt (float):                  Timestep in minutes, defaults to 0.01 minutes
    minutes (int):               Minutes the model runs for, defaults to 60 minutes
    
    Returns
    Dictionary of numpy arrays:  Returns dictionary of numpy arrays with initial concentrations 
                                 for each array
                                 
    Examples:
    
    """
    arrays = {}
    
    if len(species) != len(init_conc):
        raise ValueError('Please provide initial concentrations for each species.')
    
    for i in range(len(species)):
        arrays[species[i]] = np.zeros(int(minutes/dt))
        arrays[species[i]][i] = init_conc[i]
    
    return arrays
    
def error(*model_values, dt, minutes):
    """
    A function that calculates the relative error of the model values against data values extracted from
    figure 1 in Tsutsumi et al.
    
    Args
    model_values (numpy 2d array): Array containing model values to be assessed.
    dt (float):                    Timestep
    minutes (integer):             Number of minutes model ran
    
    """
    
    WT_model, short_model = model_values
    
    ts = np.linspace(0, minutes, int(minutes/dt))
    
    WT = np.interp(time, ts, WT_model)
    short = np.interp(0, minutes, int(minutes/dt))
    
    model = np.array([[WT],
                     [short]])
    
    return np.sum(np.power((data - model),2))