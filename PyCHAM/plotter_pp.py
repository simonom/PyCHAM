'''plots results for the particle-phase stuff'''
# simulation results are represented graphically:
# plots results for the total particle-phase concentration temporal profiles of specified components
# also the temporal profile or non-seed and non-water particle-phase

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
	# caller - marker for whether PyCHAM (0 for individual components, 3 for 
	#	total excluding seed and water, 4 for top contributors to 
	#	particle-phase, 5 for particle surface area concentration, 
	#	6 for seed surface area concentration) 
	#	or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_MW, Nwet, comp_names, y_MV, _, wall_on, space_mode, 
		_, _, _, PsatPa, OC, H2Oi, seedi, _, _, group_indx, _) = retr_out.retr_out(dir_path)
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0 or caller >= 3):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))

	if (comp_names_to_plot): # if component names specified
	
		# concentration plot ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			if (comp_names_to_plot[i].strip() == 'H2O'):
				indx_plot = [H2Oi]
				indx_plot = np.array((indx_plot))
			if (comp_names_to_plot[i].strip() == 'RO2'):
				indx_plot = (np.array((group_indx['RO2i'])))
			if (comp_names_to_plot[i].strip() == 'RO'):
				indx_plot = (np.array((group_indx['ROi'])))
				
			if (comp_names_to_plot[i].strip() != 'H2O' and comp_names_to_plot[i].strip() != 'RO2' and comp_names_to_plot[i].strip() != 'RO'):
				try: # will work if provided components were in simulation chemical scheme
					# get index of this specified component, removing any white space
					indx_plot = [comp_names.index(comp_names_to_plot[i].strip())]
					indx_plot = np.array((indx_plot))
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
			
			# particle-phase concentrations of all components (# molecules/cm3)
			if (wall_on == 1): # wall on
				ppc = yrec[:, num_comp:-num_comp]
			if (wall_on == 0): # wall off
				ppc = yrec[:, num_comp::]
			
			# particle-phase concentration of this component over all size bins (# molecules/cm3)
			conc = np.zeros((ppc.shape[0], num_sb-wall_on))
			for indxn in indx_plot: # loop through the indices
				concf = ppc[:, indxn::num_comp]
				# particle-phase concentration (ug/m3)
				conc[:, :] += ((concf/si.N_A)*y_MW[indxn])*1.e12
			
			if (comp_names_to_plot[i].strip() != 'RO2' and comp_names_to_plot[i].strip() != 'RO'): # if not a sum
			
				# plot this component
				ax0.plot(timehr, conc.sum(axis = 1), '+', linewidth = 4., label = str(str(comp_names[indx_plot[0]])+' (particle-phase)'))
			if (comp_names_to_plot[i].strip() == 'RO2'): # if is the sum of organic peroxy radicals
				ax0.plot(timehr, conc.sum(axis = 1), '-+', linewidth = 4., label = str(r'$\Sigma$RO2 (particle-phase)'))
			
			if (comp_names_to_plot[i].strip() == 'RO'): # if is the sum of organic alkoxy radicals
				ax0.plot(timehr, conc.sum(axis = 1), '-+', linewidth = 4., label = str(r'$\Sigma$RO (particle-phase)'))
				
		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		# end of gas-phase concentration sub-plot ---------------------------------------
	
	# if called by button to plot temporal profile of total particle-phase concentration 
	# excluding water and seed
	if (caller == 3):
		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		# zero water and seed
		ppc[:, H2Oi::num_comp] = 0.
		for i in seedi: # loop through seed components
			ppc[:, i::num_comp] = 0.
		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, num_sb-wall_on))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))
		# convert from # molecules/cm3 to ug/m3
		ppc = (ppc/si.N_A)*y_mwt*1.e12
		# sum over components and size bins
		ppc = np.sum(ppc, axis=1)
		
		# plot
		ax0.plot(timehr, ppc, '+', linewidth = 4., label = 'total particle-phase excluding seed and water')
		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

	# if called by button to plot top contributors to particle-phase
	if (caller == 4):
		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, (num_sb-wall_on)))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))

		# convert to ug/m3 from # molecules/cm3
		ppc = ((ppc/si.N_A)*y_mwt)*1.e12

		# sum particle-phase concentration per component over time (ug/m3)
		ppc_t = (np.sum(ppc, axis=0)).reshape(1, -1)
		
		# sum particle-phase concentrations over size bins (ug/m3), 
		# but keeping components separate
		ppc_t = np.sum(ppc_t.reshape(num_sb-wall_on, num_comp), axis=0)
		
		# convert to mass contributions
		ppc_t = ((ppc_t/np.sum(ppc_t))*100.).reshape(-1, 1)
		
		# rank in ascending order
		ppc_ts = np.flip(np.sort(np.squeeze(ppc_t)))
		
		# take just top number as supplied by user
		ppc_ts = ppc_ts[0:self.e300r]
		
		# sum concentrations over size bins and components
		ppc_sbc = np.sum(ppc, axis=1)

		# loop through to plot
		for i in range(self.e300r):
			
			# get index
			indx = np.where(ppc_t == ppc_ts[i])[0][0]
			# get name
			namei = comp_names[indx]
			# sum for this component over size bins (ug/m3)
			ppci = np.sum(ppc[:, indx::num_comp], axis=1)
			# get contribution (%)
			ppci = (ppci[ppc_sbc>0.]/ppc_sbc[ppc_sbc>0.])*100.

			ax0.plot(timehr[ppc_sbc>0.], ppci, '-+', linewidth = 4., label = namei)
	
		ax0.set_ylabel(r'Contribution to particle-phase mass concentration (%)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)
	
	# if called by button to plot total particle surface area concentration (m2/m3)
	if (caller == 5):
		# get surface area of single particles per size bin 
		# at all times (m2), note radius converted from um 
		# to m
		asp = 4.*np.pi*(x*1.e-6)**2.
		# integrate over all particles in a size bin, 
		# note conversion from /cm3 to /m3 (m2/m3)
		asp = asp*(Nwet*1.e6)	
		# sum over size bins (m2/m3)
		asp = np.sum(asp, axis=1)
		# plot
		ax0.plot(timehr, asp, '-+', linewidth = 4.)
		ax0.set_ylabel(r'Total particle-phase surface area concentration ($\rm{m^{2}\,m^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	# if called by button to plot seed particle surface area concentration (m2/m3)
	if (caller == 6):

		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on == 1): # wall on
			ppc = yrec[:, num_comp:-num_comp]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]
		
		# empty results array for seed particle volume concentrations per size bin
		ppcs = np.zeros((ppc.shape[0], (num_sb-wall_on)))

		for i in seedi: # loop through seed indices
			# get just seed component particle-phase volume concentration 
			# (cm3/cm3)
			ppcs[:, :] += (ppc[:, i::num_comp]/si.N_A)*(np.array((y_MV))[i])
		
		# convert total volume to volume per particle (cm3)
		ppcs[Nwet>0] = ppcs[Nwet>0]/Nwet[Nwet>0]
		
		# convert volume to radius (cm)
		ppcs = ((3.*ppcs)/(4.*np.pi))**(1./3.)
		
		# convert cm to m
		ppcs = ppcs*1.e-2
		# get surface area over all particles in a size bin (m2/m3), 
		# note conversion from m2/cm3 to m2/m3
		ppcs = (4.*np.pi*ppcs**2.)*Nwet*1.e6
		# sum over all size bins (m2/m3)
		ppcs = np.sum(ppcs, axis=1)
		
		# plot
		ax0.plot(timehr, ppcs, '-+', linewidth = 4.)
		ax0.set_ylabel(r'Seed particle-phase surface area concentration ($\rm{m^{2}\,m^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	# display
	if (caller == 2):
		plt.show()	

	return()