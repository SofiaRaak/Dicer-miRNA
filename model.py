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





def conc_change(theta):
    """
    
    """