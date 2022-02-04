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


def test():
    print("Module split_model.py successfully loaded.")

#split model
def conc_change(theta):    
    k1, k2, k3 = theta
    
    k1 = np.log(k1)
    k2 = np.log(k2)
    k3 = np.log(k3)
    
    #empty arrays
    WT = np.zeros(int(minutes/dt))
    short = np.zeros(int(minutes/dt))
    dicer1 = np.zeros(int(minutes/dt))
    dicer2 = np.zeros(int(minutes/dt))
    WT_dicer1 = np.zeros(int(minutes/dt))
    short_dicer2 = np.zeros(int(minutes/dt))
    mirna1 = np.zeros(int(minutes/dt))
    mirna2 = np.zeros(int(minutes/dt))
    
    #set init values
    WT[0] = WT_init
    short[0] = short_init
    dicer1[0] = dicer_init
    dicer2[0] = dicer_init
    
    for i in range(1, int(minutes/dt)):
        if WT[i-1] < 0:
            print(f'{i-1}: WT below 0: {WT[i-1]}')
        elif short[i-1] < 0:
            print(f'{i-1}: short below 0: {short[i-1]}')
        elif dicer1[i-1] < 0:
            print(f'{i-1}: dicer1 0: {dicer1[i-1]}')
        elif dicer2[i-1] < 0:
            print(f'{i-1}: dicer2 0: {dicer2[i-1]}')
        elif WT_dicer1[i-1] < 0:
            print(f'{i-1}: WT_dicer1 below 0: {WT_dicer1[i-1]}')
        elif short_dicer2[i-1] < 0:
            print(f'{i-1}: short_dice2r below 0: {short_dicer2[i-1]}')
        elif mirna1[i-1] < 0:
            print(f'{i-1}: mirna1 below 0: {mirna1[i-1]}')
        elif mirna2[i-1] < 0:
            print(f'{i-1}: mirna2 below 0: {mirna2[i-1]}')
        
        #wild type pre-mirna loop length
        WT[i] = WT[i-1] + dt*(WT_dicer1[i-1]*k_1 - WT[i-1]*dicer1[i-1]*k1)
        dicer1[i] = dicer1[i-1] + dt*(WT_dicer1[i-1]*(k_1 + k3) - dicer1[i-1] * WT[i-1] * k1)
        mirna1[i] = mirna1[i-1] + dt*(WT_dicer1[i-1]*k3)
        #mirna1[i] = (WT[0] + WT_dicer1[0]) - (WT_dicer1[i] + WT[i])
        
        #short loop pre-mirna
        short[i] = short[i-1] + dt*(short_dicer2[i-1]*k_2 - short[i-1] * dicer2[i-1]*k2)
        dicer2[i] = dicer2[i-1] + dt*(short_dicer2[i-1]*(k_2 * k3) - dicer2[i-1] * short[i-1] * k2)
        mirna2[i] = mirna2[i-1] + dt*(short_dicer2[i-1]*k3)
        #mirna2[i] = (short[0] +  short_dicer2[0]) - (short[i] - short_dicer2[i])
        
        return WT, short 


def frac_diced(theta):
    WT, short = conc_change(theta)
    
    WT_diced = np.zeros(int(minutes/dt))
    short_diced = np.zeros(int(minutes/dt))
    
    for i in range(1, int(minutes/dt)):
        WT_diced[i] = (WT[0] - WT[i]) / WT[0]
        short_diced[i] = (short[i] - short[0]) / short[0]
    
    return WT_diced, short_diced


def error(theta):
    WT_diced, short_diced = frac_diced(theta)
    
    ts = np.linspace(0, minutes, int(minutes/dt))
    
    WT_d = np.interp(time, ts, WT_diced)
    short_d = np.interp(time, ts, short_diced)
    
    model_values = np.array([[WT_d],
                            [short_d]])

    return np.sum((data_values - model_values ) ** 2)