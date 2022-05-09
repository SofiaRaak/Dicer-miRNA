import matplotlib.pyplot as plt
from scipy.optimize import minimize
import cma
import os
import utils
import model
import numpy as np
import params

def ErrorODE(theta):
    """
    Error function for ODE model based on error function described in utils module
    """
    WT_diced, short_diced, ts = model.frac_diced_ODE(model.theta)
    
    model_values =  np.array([WT_diced, short_diced])
    
    return utils.error_ODE(model_values = model_values, ts = ts)
    
    
def Error(theta):
    """
    Error function for stepwise model based on error function described in utils module
    """
    
    WT_diced, short_diced = model.frac_diced(model.theta)
    
    model_values = np.array([WT_diced, short_diced])
    
    return utils.error(model_values = model_values, dt = params.dt, minutes = params.minutes)



print(np.exp(model.theta))

print(utils.data)

print(ErrorODE(model.theta))

res1 = minimize(ErrorODE, model.theta, method = 'Nelder-Mead', tol=1e-6) #nm

print(res1)

print(r'\n\n')

print(np.exp(model.theta))

print(Error(model.theta))

res2 = minimize(Error, model.theta, method = 'Nelder-Mead', tol=1e-6)

print(res2)