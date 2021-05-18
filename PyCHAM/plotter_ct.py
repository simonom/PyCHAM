'''plots results for the change tendency temporal profiles of specified components'''
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
		y_mw, _, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, _, _, _, _, _) = retr_out.retr_out(dir_path)
	
	# loop through components to plot to check they are available
	for comp_name in (comp_names_to_plot):
		
		fname = str(dir_path+ '/' +comp_name +'_rate_of_change')
		try: # try to open
			dydt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header	
		except:
			mess = str('Please note, a change tendency record for the component ' + str(comp_name) + ' was not found, was it specified in the tracked_comp input of the model variables file?  Please see README for more information.')
			self.l203a.setText(mess)
			
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
	
	
	# if all files are available, then proceed without error message
	mess = str('')
	self.l203a.setText(mess)
			
	if (self.bd_pl < 3):
		self.l203a.setStyleSheet(0., '0px solid red', 0., 0.)
		self.bd_pl == 3
	
	
	# prepare figure
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	
	for comp_name in (comp_names_to_plot): # loop through components to plot
		
		ci = comp_names.index(comp_name) # get index of this component
		
		# note that penultimate column in dydt is gas-particle 
		# partitioning and final column is gas-wall partitioning, whilst
		# the first row contains chemical reaction numbers
		# extract the change tendency due to gas-particle partitioning
		gpp = dydt[1::, -2]
		# extract the change tendency due to gas-wall partitioning
		gwp = dydt[1::, -1]
		# sum chemical reaction gains
		crg = np.zeros((dydt.shape[0]-1, 1))
		# sum chemical reaction losses
		crl = np.zeros((dydt.shape[0]-1, 1))
		for ti in range(dydt.shape[0]-1): # loop through times
			indx = dydt[ti+1, 0:-2] > 0 # indices of reactions that produce component
			crg[ti] = dydt[ti+1, 0:-2][indx].sum()
			indx = dydt[ti+1, 0:-2] < 0 # indices of reactions that lose component
			crl[ti] = dydt[ti+1, 0:-2][indx].sum()
			 
		# convert change tendencies from molecules/cc/s to ug/m3/s
		gpp = ((gpp/si.N_A)*y_mw[ci])*1.e12
		gwp = ((gwp/si.N_A)*y_mw[ci])*1.e12
		crg = ((crg/si.N_A)*y_mw[ci])*1.e12
		crl = ((crl/si.N_A)*y_mw[ci])*1.e12
			 
		# plot temporal profiles of change tendencies due to chemical 
		# reaction production and loss, gas-particle partitioning and gas-wall partitioning
		ax0.plot(timehr, gpp, label = str('gas-particle partitioning '+ comp_name))
		ax0.plot(timehr, gwp, label = str('gas-wall partitioning '+ comp_name))
		ax0.plot(timehr, crg, label = str('chemical reaction gain '+ comp_name))
		ax0.plot(timehr, crl, label = str('chemical reaction loss '+ comp_name))
		ax0.yaxis.set_tick_params(direction = 'in')
		
		ax0.set_title('Change tendencies, where a tendency to decrease \ngas-phase concentrations is treated as negative')
		ax0.set_xlabel('Time through experiment (hours)')
		ax0.set_ylabel('Change tendency ($\mathrm{\mu g\, m^{-3}\, s^{-1}}$)')
		
		ax0.yaxis.set_tick_params(direction = 'in')
		ax0.xaxis.set_tick_params(direction = 'in')
		
		ax0.legend()
			 

	return()