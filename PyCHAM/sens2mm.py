'''solving the sensitivity (Hz/ppt) of instrument to molar mass (g/mol)'''
# module to estimate the sensitivity of an instrument to the molar mass of components, for example a Chemical Ionisiation Mass Spectrometer
# File Created at 2021-09-20 12:00:11.297818

import numpy as np

# function for sensitivity
def sens2mm(caller, y_MW):
	
	# inputs: -----------------
	# caller - flag for the calling function
	# y_MW - molar mass (g/mol) of components in question
	# ---------------------------
	
	fac_per_comp = 1. # sensitivity (Hz/ppt) per molar mass (g/mol) 
	fac_per_comp = np.array((fac_per_comp)).reshape(-1) # reshape 
	if (len(fac_per_comp) == 1): # if just a single value then tile across components 
		fac_per_comp = np.tile(fac_per_comp, len(y_MW)) # if just a single value then tile across components 
	
	if (caller == 3): # called on to plot sensitivity to molar mass
		import matplotlib.pyplot as plt 
		plt.ion()
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
		ax0.plot(y_MW, fac_per_comp)
		ax0.set_title('Sensitivity of instrument to molar mass')
		ax0.set_ylabel('Sensitivity (fraction (0-1))', size = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.set_xlabel('Molar Mass ($\mathrm{g\,mol^{-1}}$)', fontsize=14)
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	
	return(fac_per_comp)