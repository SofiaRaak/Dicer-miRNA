"""
A selection of functions required for running and optimisation of the model

"""

import numpy as np
import math




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
    
    for i in species:
        arrays[species[i]] = np.zeros(int(minutes/dt))
        arrays[species[i]][i] = inic_conc[i]
    
    return arrays
    
