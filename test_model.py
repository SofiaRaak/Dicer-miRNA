import pytest
import model
from model import conc_change
from model import frac_diced
from model import ODE_model
from model import frac_diced_ODE
from scipy.integrate import solve_ivp
import numpy as np
import math
import params
import tqdm

def test_model():
    WT, short = conc_change(model.theta)
    for i in tqdm.tqdm(range(len(WT))):
        assert type(WT[i]) is np.float64 
        assert type(short[i]) is np.float64 
        assert WT[i] >= 0 
        assert short[i] >=0
        assert WT[i] <= WT[0]
        assert short[i] <= short[0]
        

def test_fractions():
    WT, short = frac_diced(model.theta)
    for i in tqdm.tqdm(range(1, len(WT))):
        assert WT[i] >= WT[i-1]
        assert short[i] >= short[i-1]
        
def test_ODE_model():
    sol = solve_ivp(ODE_model, (0, int(params.minutes)), model.init_values, args = (model.ka, model.kb, model.kc))#, t_eval = np.linspace(0, int(params.minutes), int(params.minutes/params.dt)))#, method = 'BDF')#, args = (model.theta))
    WT, short, dicer1, dicer2, WT_dicer, short_dicer,  mirna1, mirna2 = sol.y
    print(sol.message)
    for i in tqdm.tqdm(range(len(WT))):
        assert type(WT[i]) is np.float64
        assert type(short[i]) is np.float64 
        assert WT[i] >= 0 
        assert short[i] >=0
        assert WT[i] <= WT[0]
        assert short[i] <= short[0]
        
def test_fractions_ODE():
    WT, short, time = frac_diced_ODE(model.theta)
    for i in tqdm.tqdm(range(1, len(WT))):
        assert WT[i] >= WT[i-1]
        assert short[i] >= short[i-1]