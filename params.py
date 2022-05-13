## Model setup
#timesteps
dt = 0.00001
minutes = 60

Kd_wt = 25.4 #nM, experimental values
Kd_short = 147.7 #nM, experimental values

#K_d = k_off / k_on

WT_init = 1 #nM, from experimental setup
short_init = 1 #nM
dicer_init = 5 #nM
mirna_init = 0
WT_dicer_init = 0
short_dicer_init = 0

k1 = 1
k_1 = Kd_wt * k1
k2 = 1
k_2 = Kd_short * k2
k3 = 1

#ODE optimised k1, k2, k3 (np.log) 
#theta = np.array([ 0.02817454, -0.14887925, -4.01565745])