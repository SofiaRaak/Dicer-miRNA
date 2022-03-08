import pytest
from utils import generate_arrays
import numpy as np
import math

@pytest.mark.parametrize('species, init_conc, expect', [
    (['WT', 'Mutant'], [0.1, 1], {'WT': np.array([0.1, 0, 0]), 'Mutant': np.array([1,0,0])}),
    ([1, 2], [0.5, 5], {1: np.array([0.5, 0, 0]), 2: np.array([5, 0, 0])})
])
def test_generate_arrays(species, init_conc, expect):
    test = generate_arrays(species=species, init_conc=init_conc, dt = 1, minutes = 3)
    
    assert test[species[0]][0] == expect[species[0]][0] and len(test[species[-1]]) == len(expect[species[-1]])