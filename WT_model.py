import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import cma
import os

#timesteps
dt = 0.0001
minutes = 60

Kd_wt = 25.4 #nM, experimental values

#K_d = k_off / k_on

WT_init = 1 #nM, from experimental setup
dicer_init = 5 #nM
mirna = 0

k1 = 5
k_1 = Kd_wt * k1
k3 = 5

theta = np.array([k1, k3])

#data fig 1
time = np.array([0, 5, 10, 20, 40, 60])
WT_y = np.array([0, 0.11144276160503169, 0.16566679779700877, 0.23905143587726366, 0.2954956726986665, 0.2946793863099961])

def conc_change(theta):
    k1, k3 = theta
    k1 = np.log(k1)
    k3 = np.log(k3)

    WT = np.zeros(int(minutes/dt))
    dicer = np.zeros(int(minutes/dt))
    WT_dicer = np.zeros(int(minutes/dt))
    mirna = np.zeros(int(minutes/dt))

    WT[0] = WT_init
    dicer[0] = dicer_init

    for i in range(1, int(minutes/dt)):
        WT[i] = WT[i-1] + dt*(WT_dicer[i-1]*k_1 - WT[i-1]*dicer[i-1]*k1)
        dicer[i] = dicer[i-1] + dt*(WT_dicer[i-1]*k_1 + WT_dicer[i-1]*k3 - WT[i-1]*dicer[i-1]*k1)
        WT_dicer[i] = (WT[0] + dicer[0]) - (WT[i] + dicer[i])
        mirna[i] = mirna[i-1] + dt*(WT_dicer[i-1]*k3)

    return WT

def frac_diced(theta):
    WT = conc_change(theta)
    
    WT_diced = np.zeros(int(minutes/dt))

    for i in range(1, int(minutes/dt)):
        WT_diced[i] = (WT[0] - WT[i]) / WT[0]

    return WT_diced

def error(theta):
    WT_diced = frac_diced(theta)

    ts = np.linspace(0, minutes, int(minutes/dt))
    
    WT_d = np.interp(time, ts, WT_diced)

    return np.sum((WT_y - WT_d)**2)

def test():
    print('Model module successfully imported.')