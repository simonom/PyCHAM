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
	
	# no record of change tendency for final experiment time point
	timehr = timehr[0:-1]

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
	
# plot change tendency due to individual chemical reactions
def plotter_ind(caller, dir_path, comp_names_to_plot, top_num, uc, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# top_num - top number of chemical reactions to plot
	# uc - units to use for change tendency
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
			dydt = np.loadtxt(fname, delimiter = ',', skiprows = 0) # skiprows = 0 skips first header	
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
		
		# prepare to store results of change tendency due to chemical reactions
		res = np.zeros((dydt.shape[0], dydt.shape[1]-2))
		res[:, :] = dydt[:, 0:-2] # get chemical reaction numbers and change tendencies
		
		if (uc == 0):
			Cfaca = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list
			# convert change tendencies from molecules/cc/s to ppb
			res[1::, :] = (res[1::, :]/Cfaca[1::])
			ct_units = str('(ppb/s)')
		if (uc == 1):			
			# convert change tendencies from molecules/cc/s to ug/m3/s
			res[1::, :] = ((res[1::, :]/si.N_A)*y_mw[ci])*1.e12
			ct_units = str('(' + u'\u03BC' + 'g/m' +u'\u00B3' + '/s)')
		if (uc == 2):
			# keep change tendencies as molecules/cc/s
			res[1::, :] = res[1::, :]
			ct_units = str('\n(' + u'\u0023' + ' molecules/cm' +u'\u00B3' + '/s)')
		
		# identify most active chemical reactions
		# first sum total change tendency over time (ug/m3/s)
		res_sum = np.abs(np.sum(res[1::, :], axis=0))
		 
		# sort in ascending order
		res_sort = np.sort(res_sum)
		
		if (len(res_sort) < top_num[0]):

			# if less reactions are present than the number requested inform user
			mess = str('Please note that ' + str(len(res_sort)) + ' relevant reactions were found, although a maximum of ' + str(top_num[0]) + ' were requested by the user.')
			self.l203a.setText(mess)

		for cnum in range(np.min([top_num[0], len(res_sort)])): # loop through chemical reactions
			
			# identify this chemical reaction
			cindx = np.where((res_sort[-(cnum+1)] ==  res_sum) == 1)[0]
			for indx_two in (cindx):
				# plot, note the +1 in the label to bring label into MCM index
				ax0.plot(timehr[0:-1], res[1::, indx_two], label = str(int(res[0, indx_two])+1))
		
		ax0.yaxis.set_tick_params(direction = 'in')
		
		ax0.set_title('Change tendencies, where a tendency to decrease \ngas-phase concentrations is negative')
		ax0.set_xlabel('Time through experiment (hours)')
		ax0.set_ylabel(str('Change tendency ' + ct_units))
		
		ax0.yaxis.set_tick_params(direction = 'in')
		ax0.xaxis.set_tick_params(direction = 'in')
		
		ax0.legend()
			 

	return()