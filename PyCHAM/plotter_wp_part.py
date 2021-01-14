'''plots results for the wall-phase (from particle partitioning to wall) temporal profiles of specified components'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import retr_out
import numpy as np
import scipy.constants as si

def plotter(caller, dir_path, comp_names_to_plot):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, comp_names, 
		y_MW, _, _, y_MV, _, wall_on, space_mode, _, _, yrec_p2w) = retr_out.retr_out(dir_path)
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	if (comp_names_to_plot): # if component names specified
	
		# plotting section ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			# get index of this specified component, removing any white space
			indx_plot = comp_names.index(comp_names_to_plot[i].strip())
			if (wall_on == 1):
				# total concentration on wall (from particle deposition to wall) (molecules/cc)
				conc = (yrec_p2w[:, indx_plot::num_comp]).sum(axis = 1)
				
			else:
				print('Wall not considered in this simulation')
				return()
				
			# concentration in ug/m3
			conc = ((conc/si.N_A)*y_MW[indx_plot])*1.e12
			# plot this component
			ax0.plot(timehr, conc, '+', linewidth = 4., label = str(str(comp_names[indx_plot ]+' (wall (from particle deposition to wall))')))

		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		# end of gas-phase concentration sub-plot ---------------------------------------
	

	# display
	if (caller == 2):
		plt.show()	

	return()