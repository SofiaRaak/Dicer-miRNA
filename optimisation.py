import matplotlib.pyplot as plt
from scipy.optimize import minimize
import cma
import os
import utils
import model
import numpy as np

def ErrorODE(theta):
    """
    Error function for ODE model based on error function described in utils module
    """
    WT_diced, short_diced, ts = model.frac_diced_ODE(model.theta)
    
    model_values =  np.array([WT_diced, short_diced])
    
    return utils.error_ODE(model_values = model_values, ts = ts)
    

print(model.theta)

print(ErrorODE(model.theta))

res1 = minimize(ErrorODE, model.theta)#, method = 'BFGS') #nm

print(res1)