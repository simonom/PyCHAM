'''solving the size-dependent loss rate of particles in instrument inlet'''
# module to estimate loss rate (fraction /s) of particles during passage through an instrument inlet 
# File Created at 2021-08-16 12:38:02.687302

import numpy as np

# function for loss rate
def inlet_loss_func(Dp_all):
	
	# inputs: -------------------
	# Dp_all - diameter (um) of size bins (columns) per time step (rows)
	# -----------------------------
	
	sd_lrate = np.zeros((Dp_all.shape[0], Dp_all.shape[1]))
	try: # in case loss function string is acceptable
		# estimate loss rate (fraction/s)
		for it in range(Dp_all.shape[0]):
			Dp = Dp_all[it, :]
			sd_lrate[it, :] = 0.
	except: # in case of issue
		sd_lrate = 'Error, function of loss rate of particles during passage through inlet failed, please revise'
	
	return(sd_lrate)