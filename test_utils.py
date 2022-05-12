import pytest
import utils
from utils import generate_arrays
from utils import error
from utils import error_ODE
import numpy as np
import math

minutes = 60
dt1 = 1
ts1 = np.linspace(0, minutes, int(minutes/dt1))
mock1a = np.random.random(int(minutes/dt1))
mock1a = np.sort(mock1a)
mock1a = np.interp(utils.time, ts1, mock1a)
mock1b = np.random.random(int(minutes/dt1))
mock1b = np.sort(mock1b)
mock1b = np.interp(utils.time, ts1, mock1b)
mock1 = np.array([[mock1a],
                  [mock1b]])
expect1 = np.sum(np.power((utils.data - mock1), 2))

dt2 = 0.1
ts2 = np.linspace(0, minutes, int(minutes/dt2))
mock2a = np.random.random(int(minutes/dt2))
mock2a = np.sort(mock2a)
mock2a = np.interp(utils.time, ts2, mock2a)
mock2b = np.random.random(int(minutes/dt2))
mock2b = np.sort(mock2b)
mock2b = np.interp(utils.time, ts2, mock2b)
mock2 = np.array([[mock2a],
                  [mock2b]])
expect2 = np.sum(np.power((utils.data - mock2), 2))

@pytest.mark.parametrize('species, init_conc, expect', [
    (['WT', 'Mutant'], [0.1, 1], {'WT': np.array([0.1, 0, 0]), 'Mutant': np.array([1,0,0])}),
    ([1, 2], [0.5, 5], {1: np.array([0.5, 0, 0]), 2: np.array([5, 0, 0])})
])
def test_generate_arrays(species, init_conc, expect):
    test = generate_arrays(species=species, init_conc=init_conc, dt = 1, minutes = 3)
    
    assert test[species[0]][0] == expect[species[0]][0] and len(test[species[-1]]) == len(expect[species[-1]])
    
@pytest.mark.parametrize('model_values, dt, minutes, expect', [
    (mock1, dt1, minutes, expect1),
    (mock2, dt2, minutes, expect2),
])

def test_error(model_values, dt, minutes, expect):
    assert error(model_values=model_values, dt=dt, minutes=minutes) == expect
    

@pytest.mark.parametrize('model_values, ts, expect',[
    (mock1, ts1, expect1),
    (mock2, ts2, expect2),
])

def test_ODE_error(model_values, ts, expect):
    assert error_ODE(model_values=model_values, ts=ts) == expect