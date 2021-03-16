'''solving the size-dependent loss rate of particles in instrument inlet'''
# module to estimate loss rate (fraction /s) of particles during passage through an instrument inlet 
# File Created at 2021-03-15 11:07:03.251379

import numpy as np

# function for loss rate
def inlet_loss_func(Dp_all):
	
	# inputs: -------------------
	# Dp_all - diameter (um) of size bins (columns) per time step (rows)
	# -----------------------------
	
	sd_lrate = np.zeros((Dp_all.shape[0], Dp_all.shape[1]))
	# estimate loss rate (fraction/s)
	for it in range(Dp_all.shape[0]):
		Dp = Dp_all[it, :]
		sd_lrate[it, :] = np.append(10.**(-1.5-0.5*np.log10(Dp[Dp<=1.e-1])), 10.**(-0.9+0.1*np.log10(Dp[Dp>1.e-1])))
	
	return(sd_lrate)