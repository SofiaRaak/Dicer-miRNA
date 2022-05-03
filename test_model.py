import pytest
import model
from model import conc_change
from model import frac_diced
from model import ODE_model
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
        

def test_fractions():
    WT, short = frac_diced(model.theta)
    for i in tqdm.tqdm(range(1, len(WT))):
        assert WT[i] >= WT[i-1]
        assert short[i] >= short[i-1]
        
def test_ODE_model():
    sol = solve_ivp(ODE_model, (0, int(params.minutes)), model.init_values, t_eval = np.linspace(0, int(params.minutes), int(params.minutes/params.dt)))#, args = (model.theta))
    WT, WT_dicer, dicer1, mirna1, short, short_dicer, dicer2, mirna2 = sol.y
    for i in tqdm.tqdm(range(len(WT))):
        assert type(WT[i]) is np.float64
        assert type(short[i]) is np.float64 
        assert WT[i] >= 0 
        assert short[i] >=0