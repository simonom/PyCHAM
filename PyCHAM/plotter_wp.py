'''plots results for the wall-phase (from gas-wall partitioning) temporal profiles of specified components'''
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

def plotter(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, _, _, _) = retr_out.retr_out(dir_path)
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	if (comp_names_to_plot): # if component names specified
	
		# concentration plot ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			if comp_names_to_plot[i].strip() == 'H2O':
				indx_plot = H2Oi[0]
			else:
				try: # will work if provided components were in simulation chemical scheme
					# get index of this specified component, removing any white space
					indx_plot = comp_names.index(comp_names_to_plot[i].strip())
				except:
					self.l203a.setText(str('Component ' + comp_names_to_plot[i] + ' not found in chemical scheme used for this simulation'))
					# set border around error message
					if (self.bd_pl == 1):
						self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
						self.bd_pl = 2
					else:
						self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
						self.bd_pl = 1
					
				plt.ioff() # turn off interactive mode
				plt.close() # close figure window
				return()
			
			if (wall_on == 1):
				# wall-phase concentration (from gas-particle partitioning) (molecules/cc)
				conc = yrec[:, ((num_asb+1)*num_comp)+indx_plot]
			else:
				self.l203a.setText(str('Wall not considered in this simulation'))
				
				# set border around error message
				if (self.bd_pl == 1):
					self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
					self.bd_pl = 2
				if (self.bd_pl >= 2):
					self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
					self.bd_pl = 1

				plt.ioff() # turn off interactive mode
				plt.close() # close figure window
				return()
				
			# concentration in units of ug/m3
			conc = ((conc/si.N_A)*y_MW[indx_plot])*1.e12
			# plot this component
			ax0.plot(timehr, conc, '+', linewidth = 4., label = str(str(comp_names[indx_plot ]+' (wall (from gas-wall partitioning))')))

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