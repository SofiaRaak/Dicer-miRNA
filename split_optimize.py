#import modulefinder
import split_model
#from importlib import reload
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import cma
import os

#reload(model)

split_model.test()

if sys.argv[1] == 'cma':
    try:
        res = cma.fmin(split_model.error(split_model.theta), split_model.theta, 2)
    except TypeError:
        os.system('rename outcmaes outcmaes_del')
        res = cma.fmin(split_model.error(split_model.theta), split_model.theta, 2)

if sys.argv[1] == 'nm':
    res = minimize(split_model.error, split_model.theta, method = 'nelder-mead')

if sys.argv[1] == 'bfgs':
    res = minimize(split_model.error, split_model.theta, method = 'bfgs')

if sys.argv[1] == 'powell':
    res = minimize(split_model.error, split_model.theta, method = 'powell')

else:
    print("Input not recognised! Valid inputs: 'cma', 'nm', 'bfgs', 'powell'.")

print(res)