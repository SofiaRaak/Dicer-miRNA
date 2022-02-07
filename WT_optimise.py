import WT_model as model
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import cma
import os

model.test()

if sys.argv[1] == 'cma':
    try:
        res = cma.fmin(model.error(model.theta), model.theta, 2)
    except TypeError:
        os.system('rename outcmaes outcmaes_del')
        res = cma.fmin(model.error(model.theta), model.theta, 2)

if sys.argv[1] == 'nm':
    res = minimize(model.error, model.theta, method = 'nelder-mead')

if sys.argv[1] == 'bfgs':
    res = minimize(model.error, model.theta, method = 'bfgs')

if sys.argv[1] == 'powell':
    res = minimize(model.error, model.theta, method = 'powell')

else:
    print("Input not recognised! Valid inputs: 'cma', 'nm', 'bfgs', 'powell'.")

print(res)