import pytest
import model
from model import conc_change
from model import frac_diced
from model import ODE_model
from scipy.integrate import solve_ivp
import numpy as np
import math
import params

def test_model():
    model_values = conc_change(model.theta)
    for i in range(len(model_values[0])):
        assert type(model_values[0][i]) is float and type(model_values[1]) is float and model_values[0][i] >= 0 and model_values[0][1] >=0
        

def test_fractions():
    fracs = frac_diced(model.theta)
    for i in range(1, len(fracs[0])):
        assert fracs[0][i] >= fracs[0][i-1] and fracs[1][i] >= fracs[1][i-1]
        
def test_ODE_model():
    sol = solve_ivp(ODE_model, (0, int(params.minutes)), model.init_values, t_eval = np.linspace(0, int(params.minutes), int(params.minutes/params.dt)))
    model_values = sol.y
    for i in range(len(model_values[0])):
        assert type(model_values[0][i]) is float and type(model_values[1]) is float and model_values[0][i] >= 0 and model_values[0][1] >=0