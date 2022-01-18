## setup environment

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.optimize import curve_fit
import cma
import os

##setup initial values for model
#timesteps
dt = 0.0001
minutes = 60

Kd_wt = 25.4 #nM, experimental values
Kd_short = 147.7 #nM, experimental values

#K_d = k_off / k_on

WT_init = 1 #nM, from experimental setup
short_init = 1 #nM
dicer_init = 5 #nM
mirna = 0

k1 = 5
k_1 = Kd_wt * k1
k2 = 5
k_2 = Kd_short * k2
k3 = 5

theta = [k1, k2, k3]

#data fig 1
time = np.array([0, 5, 10, 20, 40, 60])
WT_y = np.array([0, 0.11144276160503169, 0.16566679779700877, 0.23905143587726366, 0.2954956726986665, 0.2946793863099961])
short_y = np.array([0, 0.0033684107002276975, 0.007599822974028003, 0.010019177812737812, 0.009603658536577298, 0.01242378048779691])

data_values = np.array([[WT_y],
                       [short_y]])

##create tester function to ensure model module loaded
def test():
    print("Module model.py successfully loaded.")

##model
def conc_change(theta):    
    k1, k2, k3 = theta
    k1 = k1**2
    k2 = k2**2
    k3 = k3**2
    
    #empty arrays
    WT = np.zeros(int(minutes/dt))
    short = np.zeros(int(minutes/dt))
    dicer = np.zeros(int(minutes/dt))
    WT_dicer = np.zeros(int(minutes/dt))
    short_dicer = np.zeros(int(minutes/dt))
    mirna = np.zeros(int(minutes/dt))
    
    #set init values
    WT[0] = WT_init
    short[0] = short_init
    dicer[0] = dicer_init
    
    for i in range(1, int(minutes/dt)):
        if WT[i-1] < 0:
            print(f'{i-1}: WT below 0: {WT[i-1]}')
        elif short[i-1] < 0:
            print(f'{i-1}: short below 0: {short[i-1]}')
        elif dicer[i-1] < 0:
            print(f'{i-1}: dicer 0: {dicer[i-1]}')
        elif WT_dicer[i-1] < 0:
            print(f'{i-1}: WT_dicer below 0: {WT_dicer[i-1]}')
        elif short_dicer[i-1] < 0:
            print(f'{i-1}: short_dicer below 0: {short_dicer[i-1]}')
        elif mirna[i-1] < 0:
            print(f'{i-1}: mirna below 0: {mirna[i-1]}')
            
        WT[i] = WT[i-1] + dt * (WT_dicer[i-1] * k_1 - WT[i-1] * dicer[i-1] * k1)
        short[i] = short[i-1] + dt * (short_dicer[i-1] * k_2 - short[i-1] * dicer[i-1] * k2)
        dicer[i] = dicer[i-1] + dt * (WT_dicer[i-1] * (k_1 + k3) + short_dicer[i-1] * (k_2 + k3) - dicer[i-1] * (WT[i-1] * k1 + short[i-1] * k2))
        WT_dicer[i] = WT_dicer[i-1] + dt * (WT[i-1] * dicer[i-1] * k1 - WT_dicer[i-1] * (k_1 + k3))
        short_dicer[i] = short_dicer[i-1] + dt *(short[i-1] * dicer[i-1] * k2 - short_dicer[i-1] * (k_2 + k3))
        mirna[i] = mirna[i-1] + dt * (k3 * (WT_dicer[i-1] + short_dicer[i-1]))
    
    return WT, short, dicer, WT_dicer, short_dicer, mirna

##Fraction diced
def frac_diced(theta):
    
    WT, short, dicer, WT_dicer, short_dicer, mirna = conc_change(theta)
    
    WT_diced = np.zeros(int(minutes/dt))
    short_diced = np.zeros(int(minutes/dt))
    
    for i in range(1, int(minutes/dt)):
        WT_diced[i] = (WT[0] - WT[i]) / WT[0]
        short_diced[i] = (short[i] - short[0]) / short[0]
    
    return WT_diced, short_diced

##Error function
#create error function
def error(theta):
    WT_diced, short_diced = frac_diced(theta)
    
    ts = np.linspace(0, minutes, int(minutes/dt))
    
    WT_d = np.interp(time, ts, WT_diced)
    short_d = np.interp(time, ts, short_diced)
    
    model_values = np.array([[WT_d],
                            [short_d]])
    
    return np.sum((data_values - model_values ) ** 2)