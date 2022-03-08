import pytest
from arrays import generate_arrays
import numpy as np
import math

@pytest.mark.parameterize('species', 'init_conc', 'expect', [
    (['WT', 'Mutant'], [0.1, 1], dict('WT': array([0.1, 0, 0]),
                                     'mutant': array([1,0,0]))),
    ([1, 2], [0.5, 5], {1: array([0.5, 0, 0]), 2: array([5, 0, 0])})
])
def test_generate_arrays(species, init_conc, expect):
    assert add_arrays(species=species, init_conc=init_conc, dt = 1, minutes = 3) == expect