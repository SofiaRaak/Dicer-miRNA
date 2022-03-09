import matplotlib.pyplot as plt
from scipy.optimize import minimize
import cma
import os
import utils
import model

model_values = model.frac_diced(model.theta)

print(model.theta)

print(utils.error(model_values = model_values, dt = model.dt, minutes = model.minutes))

#res1 = minimize(utils.error(model_values = model_values, dt = model.dt, minutes = model.minutes), model.theta) #nm