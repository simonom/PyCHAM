##########################################################################################
#                                                                                        #
#    Copyright (C) 2018-2025 Simon O'Meara : simon.omeara@manchester.ac.uk               #
#                                                                                        #
#    All Rights Reserved.                                                                #
#    This file is part of PyCHAM                                                         #
#                                                                                        #
#    PyCHAM is free software: you can redistribute it and/or modify it under             #
#    the terms of the GNU General Public License as published by the Free Software       #
#    Foundation, either version 3 of the License, or (at your option) any later          #
#    version.                                                                            #
#                                                                                        #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT               #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              #
#    details.                                                                            #
#                                                                                        #
#    You should have received a copy of the GNU General Public License along with        #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                #
#                                                                                        #
##########################################################################################
'''plots results for the particle-phase stuff'''
# simulation results are represented graphically:
# plots results for the total particle-phase concentration temporal profiles of 
# specified components
# also the temporal profile or non-seed and non-water particle-phase

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import os
import numpy as np
import scipy.constants as si # for scientific constants

def plotter(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0 for individual components, 3 for 
	#	total excluding seed and water, 4 for top contributors to 
	#	particle-phase, 5 for particle surface area concentration, 
	#	6 for seed surface area concentration, 7 for top contributors
	#	to particle-phase excluding seed and water) 
	#	or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = self.ro_obj.comp_MW
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	group_indx = self.ro_obj.gi
	y_MV = np.array((self.ro_obj.comp_MV)) # cm3/mol
	
	# number of actual particle size bins
	num_asb = (self.ro_obj.nsb-self.ro_obj.wf)

	# prepare to store summed results
	if hasattr(self, 'sum_ornot_flag'):
		if (self.sum_ornot_flag == 1):
			sum_conc = np.zeros((len(timehr)))
			sum_comp_names_to_plot = comp_names_to_plot[0]

	if (caller == 0 or caller >= 3):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))

	if (comp_names_to_plot): # if component names specified
	
		ip_fail = 0 # start by assuming all requested components available
		group_flag = 0 # start by assuming single components wanted

		# concentration plot ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):

			# begin array containing indices of components to consider
			indx_plot = (np.arange(num_comp)).astype('int')
			
			if (comp_names_to_plot[i].strip() == 'H2O'):
				indx_plot = [H2Oi]
				indx_plot = np.array((indx_plot))

			if (comp_names_to_plot[i].strip() == 'RO2'):
				indx_plot = (np.array((group_indx['RO2i'])))
				group_flag = 1

			if (comp_names_to_plot[i].strip() == 'RO'):
				indx_plot = (np.array((group_indx['ROi'])))
				group_flag = 1	

			if (comp_names_to_plot[i].strip() == 'HOM'):
				indx_plot = (np.array((group_indx['HOMs'])))			
				group_flag = 1	
				
			if ('-OOH' in comp_names_to_plot[i].strip()):
				indx_plot = (np.array((group_indx['OOH'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
				
			if ('-OH' in comp_names_to_plot[i].strip()):
				indx_plot = (np.array((group_indx['OH'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1			
			
			if ('-carbonyl' in comp_names_to_plot[i].strip()):	
				indx_plot = (np.array((group_indx['carbonyl'])))			
					
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			
			if ('-NO3' in comp_names_to_plot[i].strip()):
				indx_plot = (np.array((group_indx['NO3'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			
			if ('ROOR' in comp_names_to_plot[i].strip()):
				indx_plot = (np.array((group_indx['ROOR'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1

			if ('seed' in comp_names_to_plot[i].strip()):
				indx_plot = seedi
				group_flag = 1

			# if not a group
			# check on whether this could be an individual component
			# will work if provided components were in simulation chemical scheme
			if (group_flag == 0):
				try:
					# get index of this specified component, removing any 
					# white space
					indx_plot = [comp_names.index(
						comp_names_to_plot[i].strip())]
					indx_plot = np.array((indx_plot))
					group_flag = -1
				except:
					# if not a chemical scheme component, then mark for
					# failure to find component, which will be 
					# changed if the carbon or oxygen number tests
					# below succeed
					ip_fail = 1

			if (group_flag == 0):
				if ('C' in comp_names_to_plot[i].strip() or 
					'c' in comp_names_to_plot[i].strip()):
					try:
						# index of C letter
						cindx = comp_names_to_plot[i].strip().index('C') 
					except:
						# index of c letter
						cindx = comp_names_to_plot[i].strip().index('c') 
				
					# if a number given after carbon atom
					if comp_names_to_plot[i].strip()[cindx+1].isnumeric():
						numb_indx = cindx+1 # number index finishes
						if (len(comp_names_to_plot[i].strip()) >= 
							cindx+2):
							for str_i in range(cindx+2, 
							len(comp_names_to_plot[i].strip())):
								if (comp_names_to_plot[
									i].strip()[cindx+2
									:str_i+1].isnumeric()):
									numb_indx += 1
					
						# get carbon number
						Cn = float(comp_names_to_plot[
							i].strip()[cindx+1:numb_indx+1])
					
						# get all components with this many carbons
						indx_plot_nw = np.zeros((len(indx_plot)))
						indx_plot_nw[:] = indx_plot[:]
						for indxi in indx_plot:
							# get the relevant SMILES
							SMIi = rel_SMILES[indxi]
							if ((SMIi.count('C') + 
								SMIi.count('c')) != Cn):
								indx_to_remove = np.where(
								indx_plot_nw == indxi)[0][0]
								indx_plot_nw = np.concatenate((
								indx_plot_nw[0:indx_to_remove],
								indx_plot_nw[
								indx_to_remove+1::]))
						indx_plot = (np.array((
							indx_plot_nw))).astype('int')
							
						group_flag = 1
						ip_fail = 0
						if (indx_plot.shape[0] == 0):
							ip_fail = 1

				if ('O' in comp_names_to_plot[i].strip() or 
					'o' in comp_names_to_plot[i].strip()):
					try:
						# index of O letter
						cindx = comp_names_to_plot[i].strip().index('O') 
					except:
						# index of o letter
						cindx = comp_names_to_plot[i].strip().index('o') 
				
					# if a number given after the atom letter
					if (comp_names_to_plot[i].strip()[cindx+1].isnumeric()):
						numb_indx = cindx+1 # number index finishes
						if (len(comp_names_to_plot[i].strip()) >= 
							cindx+2):
							for str_i in range(cindx+2, 
								len(comp_names_to_plot[
								i].strip())):
								if (comp_names_to_plot[
									i].strip()[
									cindx+2:str_i+1
									].isnumeric()):
									numb_indx += 1
					
						# get number of atoms of interest
						Cn = float(comp_names_to_plot[i].strip()[
							cindx+1:numb_indx+1])
					
						# get all components with this many atoms 
						# of interest
						indx_plot_nw = np.zeros((len(indx_plot)))
						indx_plot_nw[:] = indx_plot[:]
						for indxi in indx_plot:
							# get the relevant SMILES
							SMIi = rel_SMILES[indxi]
							if ((SMIi.count('O') + SMIi.count('o')) != Cn):
								indx_to_remove = np.where(indx_plot_nw == indxi)[0][0]
								indx_plot_nw = np.concatenate((indx_plot_nw[0:indx_to_remove], indx_plot_nw[indx_to_remove+1::]))
						indx_plot = (np.array((indx_plot_nw))).astype('int')
							
						ip_fail = 0
						group_flag = 1
						if (indx_plot.shape[0] == 0):
							ip_fail = 1

			
				if (ip_fail == 1):	
					# in case code thought this was a group but 
					# couldn't make the group work
					# will work if provided components were in 
					# simulation chemical scheme
					try:
						# get index of this specified component, 
						#removing any white space
						indx_plot = [comp_names.index(
							comp_names_to_plot[i].strip())]
						indx_plot = np.array((indx_plot))
						group_flag = 0
					except:

						self.l203a.setText(str('Component ' + 
						comp_names_to_plot[i] + 
						' not found in chemical scheme used ' +	
						'for this simulation'))
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
			
			if (comp_names_to_plot[i].strip() != 'H2O' and group_flag != 1):
				# will work if provided components were in simulation 
				# chemical scheme
				try:
					# get index of this specified component, 
					#removing any white space
					indx_plot = [comp_names.index(comp_names_to_plot[
						i].strip())]
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
			
			import scipy.constants as si # for scientific constants
			
			yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
			yrec[:, :] = self.ro_obj.yrec[:, :]
			
			# particle-phase concentrations of all components (# molecules/cm^3)
			if (self.ro_obj.wf > 0): # wall on
				ppc = yrec[:, self.ro_obj.nc:-self.ro_obj.nc*self.ro_obj.wf]
			if (self.ro_obj.wf == 0): # wall off
				ppc = yrec[:, self.ro_obj.nc::]
			
			# particle-phase concentration of this component over all 
			# size bins (# molecules/cm^3)
			conc = np.zeros((ppc.shape[0], self.ro_obj.nsb-self.ro_obj.wf))
			for indxn in indx_plot: # loop through the indices
				concf = ppc[:, indxn::self.ro_obj.nc]
				# particle-phase concentration (ug/m^3)
				conc[:, :] += ((concf/si.N_A)*self.ro_obj.comp_MW[indxn])*1.e12
			
			# if summing values
			if (self.sum_ornot_flag == 1):
				sum_conc += conc.sum(axis = 1)
				if (i > 0):
					sum_comp_names_to_plot = str(sum_comp_names_to_plot + ', ' + comp_names_to_plot[i].strip())
				if (i == len(comp_names_to_plot)-1):
					ax0.plot(timehr, sum_conc, '-+', linewidth = 4., label = str(r'$\Sigma$' + sum_comp_names_to_plot + ' (particle-phase)'))

			# if showing inidividual values
			if (self.sum_ornot_flag == 0):
				if (group_flag != 1): # if not a sum over a group of components
					# plot this component
					ax0.plot(timehr, conc.sum(axis = 1), '+', linewidth = 4., label = str(str(self.ro_obj.names_of_comp[indx_plot[0]]) + ' (particle-phase)'))
			
				else: # if a sum over a group of components
					ax0.plot(timehr, conc.sum(axis = 1), '-+', linewidth = 4., label = str(r'$\Sigma$' + comp_names_to_plot[i].strip() + ' (particle-phase)'))
			
			

		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		# end of particle-phase concentration sub-plot ----------------------------
	
	# if called by button to plot temporal profile of total particle-phase concentration 
	# excluding water and seed
	if (caller == 3):
		import scipy.constants as si
		# particle-phase concentrations of all components (# molecules/cm3)
		if (self.ro_obj.wf > 0): # wall on
			ppc = yrec[:, num_comp:-num_comp*self.ro_obj.wf]
		if (self.ro_obj.wf == 0): # wall off
			ppc = yrec[:, num_comp::]
		
		# zero water and seed
		ppc[:, self.ro_obj.H2O_ind::num_comp] = 0.
		for i in self.ro_obj.seed_ind: # loop through seed components
			ppc[:, i::num_comp] = 0.
		
		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((self.ro_obj.comp_MW)).reshape(1, -1), 
			(1, self.ro_obj.nsb-self.ro_obj.wf))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))
		# convert from # molecules/cm3 to ug/m3
		ppc = (ppc/si.N_A)*y_mwt*1.e12
		
		# sum over components and size bins
		ppc = np.sum(ppc, axis=1)
		
		# plot
		ax0.plot(self.ro_obj.thr, ppc, '+', linewidth = 4., 
			label = 'total particle-phase excluding seed and water')
		ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

	# if called by button to plot top contributors to particle-phase
	if (caller == 4):
		import scipy.constants as si
		# particle-phase concentrations of all components (# molecules/cm3)
		if (self.ro_obj.wf > 0): # wall on
			ppc = yrec[:, self.ro_obj.nc:-self.ro_obj.nc*self.ro_obj.wf]
		if (self.ro_obj.wf == 0): # wall off
			ppc = yrec[:, self.ro_obj.nc::]

		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((self.ro_obj.comp_MW)).reshape(1, -1), 
			(1, (self.ro_obj.nsb-self.ro_obj.wf)))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))

		# convert to ug/m3 from # molecules/cm3
		ppc = ((ppc/si.N_A)*y_mwt)*1.e12

		# sum particle-phase concentration per component over time (ug/m3)
		ppc_t = (np.sum(ppc, axis=0)).reshape(1, -1)
		
		# sum particle-phase concentrations over size bins (ug/m3), 
		# but keeping components separate
		ppc_t = np.sum(ppc_t.reshape(self.ro_obj.nsb-self.ro_obj.wf, 
			self.ro_obj.nc), axis=0)
		
		# convert to mass contributions
		ppc_t = ((ppc_t/np.sum(ppc_t))*100.).reshape(-1, 1)
		
		# rank in ascending order
		ppc_ts = np.flip(np.sort(np.squeeze(ppc_t)))
		
		# take just top number as supplied by user
		ppc_ts = ppc_ts[0:self.e300r]
		
		# sum concentrations over size bins and components
		ppc_sbc = np.sum(ppc, axis=1)

		# loop through to plot
		for i in range(np.min([self.e300r, len(self.ro_obj.names_of_comp)])):
			
			# get index
			indx = np.where(ppc_t == ppc_ts[i])[0][0]
			# get name
			namei = self.ro_obj.names_of_comp[indx]
			# sum for this component over size bins (ug/m3)
			ppci = np.sum(ppc[:, indx::self.ro_obj.nc], axis=1)
			# get contribution (%)
			ppci = (ppci[ppc_sbc>0.]/ppc_sbc[ppc_sbc>0.])*100.

			ax0.plot(self.ro_obj.thr[ppc_sbc>0.], ppci, '-+', linewidth = 4., 
				label = namei)
	
		ax0.set_ylabel(r'Contribution to particle-phase mass concentration (%)', 
			fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)
	
	# if called by button to plot top contributors to particle-phase excluding 
	# seed and water
	if (caller == 7):

		import scipy.constants as si		

		# particle-phase concentrations of all components (# molecules/cm3)
		if (self.ro_obj.wf > 0): # wall on
			ppc = yrec[:, self.ro_obj.nc:-self.ro_obj.nc*self.ro_obj.wf]
		if (self.ro_obj.wf == 0): # wall off
			ppc = yrec[:, self.ro_obj.nc::]

		for i in seedi: # loop through seed components
			ppc[:, i::self.ro_obj.nc] = 0. # zero seed

		ppc[:, self.ro_obj.H2O_ind::self.ro_obj.nc] = 0. # zero water

		# tile molar weights over size bins and times
		y_mwt = np.tile(np.array((self.ro_obj.comp_MW)).reshape(1, -1), 
			(1, (self.ro_obj.nsb-self.ro_obj.wf)))
		y_mwt = np.tile(y_mwt, (ppc.shape[0], 1))

		# convert to ug/m3 from # molecules/cm3
		ppc = ((ppc/si.N_A)*y_mwt)*1.e12

		# sum particle-phase concentration per component over time (ug/m3)
		ppc_t = (np.sum(ppc, axis=0)).reshape(1, -1)
		
		# sum particle-phase concentrations over size bins (ug/m3), 
		# but keeping components separate
		ppc_t = np.sum(ppc_t.reshape(self.ro_obj.nsb-self.ro_obj.wf, num_comp), axis=0)
		
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
			namei = self.ro_obj.names_of_comp[indx]
			# sum for this component over size bins (ug/m3)
			ppci = np.sum(ppc[:, indx::self.ro_obj.nc], axis=1)
			# get contribution (%)
			ppci = (ppci[ppc_sbc>0.]/ppc_sbc[ppc_sbc>0.])*100.

			ax0.plot(self.ro_obj.thr[ppc_sbc>0.], ppci, '-+', 
				linewidth = 4., label = namei)
	
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
		asp = 4.*np.pi*(self.ro_obj.cen_size*1.e-6)**2.
		# integrate over all particles in a size bin, 
		# note conversion from /cm3 to /m3 (m2/m3)
		asp = asp*(self.ro_obj.Nrec_wet*1.e6)	
		# sum over size bins (m2/m3)
		asp = np.sum(asp, axis=1)
		# plot
		ax0.plot(self.ro_obj.thr, asp, '-+', linewidth = 4.)
		ax0.set_ylabel(str(r'Total particle-phase surface area concentration ' +
		r'($\rm{m^{2}\,m^{-3}}$)'), fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	# if called by button to plot seed particle surface area concentration (m2/m3)
	if (caller == 6):

		import scipy.constants as si

		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on > 0): # wall on
			ppc = yrec[:, num_comp:-num_comp*wall_on]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]
		
		# empty results array for seed particle volume concentrations per size bin
		ppcs = np.zeros((ppc.shape[0], (num_sb-wall_on)))

		for i in seedi: # loop through seed indices

			# get just seed component particle-phase volume concentration 
			# (cm3/cm3)
			ppcs[:, :] += (ppc[:, i::num_comp]/si.N_A)*y_MV[i]
		
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

	# if called by button to plot generational contribution to particle-phase mass
	if (caller == 8):

		import sch_interr # for interpeting chemical scheme
		import re # for parsing chemical scheme
		import scipy.constants as si

		gen_num = [] # empty results list
		ci = 0 # component count
		sch_name = self.ro_obj.sp
		inname = self.ro_obj.vp

		f_open_eqn = open(sch_name, mode='r') # open the chemical scheme file
		# read the file and store everything into a list
		total_list_eqn = f_open_eqn.readlines()
		f_open_eqn.close() # close file

		inputs = open(inname, mode= 'r' ) # open model variables file
		in_list = inputs.readlines() # read file and store everything into a list
		inputs.close() # close file
		for i in range(len(in_list)): # loop through supplied model variables to interpret

			# ----------------------------------------------------
			# if commented out continue to next line
			if (in_list[i][0] == '#'):
				continue
			key, value = in_list[i].split('=') # split values from keys
			# model variable name - a string with bounding white space removed
			key = key.strip()
			# ----------------------------------------------------
			# formatting for chemical scheme
			if key == 'chem_scheme_markers' and (value.strip()):
				self.chem_sch_mrk = [str(i).strip() for i in (value.split(','))]

		# interrogate scheme to list equations
		self = sch_interr.sch_interr(total_list_eqn, self)	
		
		# loop through components to identify their generation
		for compi in comp_names[0:-2]: # don't include core and water
			
			comp_fin = 0 # flag for when component generation number found
			gen_num.append(0) # assume zero generation by default
			
			# if an organic molecule
			if ('C' in rel_SMILES[ci] or 'c' in rel_SMILES[ci]):
				
				# loop through reactions to identify where this first occurs,
				# assuming that if first occurrence is as a reactant 
				# it is zero generation and if as a product is >zero generation
				for eqn_step in range(len(self.eqn_list)):

					line = self.eqn_list[eqn_step] # extract this line
					# work out whether equation or reaction rate coefficient part comes first
					eqn_start = str('.*\\' +  self.chem_sch_mrk[10])
					rrc_start = str('.*\\' +  self.chem_sch_mrk[9])
					# get index of these markers, note span is the property of the match object that
					# gives the location of the marker
					eqn_start_indx = (re.match(eqn_start, line)).span()[1]
					rrc_start_indx = (re.match(rrc_start, line)).span()[1]
		
					if (eqn_start_indx > rrc_start_indx):
						eqn_sec = 1 # equation is second part
					else:
						eqn_sec = 0 # equation is first part
		
					# split the line into 2 parts: equation and rate coefficient
					# . means match with anything except a new line character., when followed by a * 
					# means match zero or more times (so now we match with all characters in the line
					# except for new line characters, so final part is stating the character(s) we 
					# are specifically looking for, \\ ensures the marker is recognised
					if eqn_sec == 1:
						eqn_markers = str('\\' +  self.chem_sch_mrk[10]+ '.*\\' +  self.chem_sch_mrk[11])
					else: # end of equation part is start of reaction rate coefficient part
						eqn_markers = str('\\' +  self.chem_sch_mrk[10]+ '.*\\' +  self.chem_sch_mrk[9])

					# extract the equation as a string ([0] extracts the equation section and 
					# [1:-1] removes the bounding markers)
					eqn = re.findall(eqn_markers, line)[0][1:-1].strip()
					
					eqn_split = eqn.split()
					eqmark_pos = eqn_split.index('=')
					# reactants with stoichiometry number and omit any photon
					reactants = [i for i in eqn_split[:eqmark_pos] if i != '+' and i != 'hv']
					# products with stoichiometry number
					products = [t for t in eqn_split[eqmark_pos+1:] if t != '+']
					
					for ri in reactants:
						# note that no spaces or other punctuation included around 
						# the component name as extracted from the equation
						if compi in ri and len(compi) == len(ri):
							# assuming that components appearing first as a reactant 
							# must be zero generation
							gen_num[ci] = 0
							# finished with this component, break out of reactant loop
							comp_fin = 1
							break

					for pi in products:
						# note that no spaces or other punctuation included around 
						# the component name as extracted from the equation
						if compi in pi and len(compi) == len(pi):
							# check which reactant has earliest generation
							rcheck = []
							for ri in reactants:
								# get its index
								rindx = comp_names.index(ri)
								# if an organic, note if inorgnic then ignore
								if ('C' in rel_SMILES[rindx] or 'c' in rel_SMILES[rindx]):
									if (rindx < ci): # if generation number of recatant already known
										rcheck.append(gen_num[rindx])
							# identify earliest generation reactant
							egen = np.min(rcheck)
							# check whether this species has a charge, i.e. is open shell, note that according 
							# to DAYLIGHT, the square brackets in a SMILES string deontes when an atom has a 
							# valence other than normal, i.e. is in an excited (radical) state: 
							# https://daylight.com/dayhtml/doc/theory/theory.smiles.html
							if ('[' in rel_SMILES[ci] or ']' in rel_SMILES[ci]): # if open-shell
								gen_num[ci] = egen
							else: # if closed-shell
								gen_num[ci] = egen+1
							# finished with this component, break out of reactant loop
							comp_fin = 1
							break

					# finished with this component, break out of equation loop, 
					# move onto next component
					if (comp_fin == 1):
						break
				

			ci += 1 # component count
		# convert list to numpy array
		gen_num = np.array(gen_num)
		# append two zeros for core and water
		gen_num = np.concatenate((gen_num, np.zeros(2)), axis=0)
		# in case we want to see generation number per component -------------
		#ax0.plot(np.arange(len(comp_names)), gen_num, '+')
		#ax0.set_xticks(np.arange(len(comp_names)))
		#ax0.set_xticklabels(comp_names, rotation = 90)
		# --------------------------------------------------------------------
		
		# empty results matrix for contribution (%) to particle-phase mass per 
		# generation per time step
		res = np.zeros((len(timehr), int(np.max(gen_num))))
		# particle-phase concentrations of all components (# molecules/cm3)
		if (wall_on > 0): # wall on
			ppc = yrec[:, num_comp:-num_comp*wall_on]
		if (wall_on == 0): # wall off
			ppc = yrec[:, num_comp::]

		y_MW = np.tile(y_MW, (len(timehr), num_sb-wall_on))

		# convert particle-phase concentrations into mass 
		# concentration (ug/m^3) from # molecules/cm^3		
		ppc[:, :] += ((ppc/si.N_A)*y_MW)*1.e12
		# zero water and seed contribution
		ppc[:, int(H2Oi)::int(num_comp)] = 0.
		for seei in seedi:
			ppc[:, int(seei)::int(num_comp)] = 0.
		# sum concentrations over all size bins 
		# and components per time step (ug/m^3)
		ppc_sum = np.sum(ppc, axis= 1)
		# non-zero indices of ppc_sum
		nzi = np.squeeze(ppc_sum != 0.)

		# tile generation number over size bins
		gen_num = np.tile(gen_num, (num_sb-wall_on))
		# loop through generations
		for gi in range(int(np.max(gen_num))):
			# get percentage contribution from this generation
			res[nzi, gi] = (np.sum(ppc[:, gen_num == gi], axis = 1)[nzi]/ppc_sum[nzi])*100.
			ax0.plot(timehr, res[:, gi], label = str('Generation ' + str(gi)))
		ax0.set_ylabel(r'Contribution to organic particle-phase mass concentration ($\%$)', fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)

		return()

	# display
	if (caller == 2):
		plt.show()	

	return()

# function for when calling from button for plotting size-segregated mass
# concentration without water
def part_mass_vs_time_sizeseg(self):

	import scipy.constants as si
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], 
		self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], 
		self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], 
		self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = self.ro_obj.comp_MW
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], 
		self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]

	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))
	
	try: # in case user-supplied values given:
		psb_dub = [float(i) for i in 
			self.e303p.text().split(',')]
	except:
		# default
		# particle size bins (upper bounds) diameter (um), 
		# ready for indexing the concentrations that will be summed
		psb_dub = [1.e-1, 5.e-1, 1.e0, 2.5e0, 1.e1, 2.e1]
	
	# convert diameters to nm from um
	for i in range(len(psb_dub)):
		psb_dub[i] = psb_dub[i]*1.e3

	# convert recorded radius to nm from um and from radius to diameter
	dbou_rec = self.ro_obj.rad*1.e3*2.
	
	# index array preparation
	psb_ind = np.zeros((len(self.ro_obj.thr), num_asb))

	# get the index of particle size bins
	for psb_dubi in psb_dub:
		
		psb_ind += dbou_rec[:, 1::]<=psb_dubi
		
	# then use this index to find the mass concentrations 
	# inside each size bound

	# concentrations of components in particles (# molecules/cm3)
	yp = yrec[:, self.ro_obj.nc:self.ro_obj.nc*(num_asb+1)]
	# zero water
	yp[:, self.ro_obj.H2O_ind::self.ro_obj.nc] = 0

	# mass concentrations of each component in each size bin (ug/m^3)
	mc = (yp/si.N_A)*np.tile(np.array((self.ro_obj.comp_MW)).reshape(1, -1), (1, num_asb))*1e12

	# mass concentrations summed across components in each size bin (ug/m^3)
	mc_psb = np.zeros((len(self.ro_obj.thr), num_asb))
	
	for i in range(num_asb): # loop through size bins
		mc_psb[:, i] = np.sum(mc[:, i*self.ro_obj.nc:(i+1)*self.ro_obj.nc], axis=1)

	# average mass per particle in this size bin (ug/m^3/particle), note that components have been
	# allocated to size bins based on their wet (with water) sizes, and we want to display
	# the mass based on dry sizes
	mc_psb[self.ro_obj.Nrec_dry>0.] = mc_psb[self.ro_obj.Nrec_dry>0.]/self.ro_obj.Nrec_dry[self.ro_obj.Nrec_dry>0.]
	
	# now get mass integrated over all particles based on the dry (no water) number concentration
	# this will ensure the mass is allocated to the correct size bin
	mc_psb = mc_psb*self.ro_obj.Nrec_wet

	# results array preparation
	psb_res = np.zeros((len(self.ro_obj.thr), len(psb_dub)))

	# temporary holder array
	mc_psb_temp = np.zeros((len(self.ro_obj.thr), num_asb))

	# loop through size categories in ascending order
	for psb_dubi in range(len(psb_dub)):
		
		# temporary holder array
		mc_psb_temp[:, :] = mc_psb[:, :]

		# zero the unwanted size bins
		mc_psb_temp[psb_ind<(len(psb_dub)-psb_dubi)] = 0.

		# get the concentrations of components in this size category (ug/m^3)
		psb_res[:, psb_dubi] = np.sum(mc_psb_temp, axis=1)
	
		ax0.plot(self.ro_obj.thr, psb_res[:, psb_dubi], label = str(str(r'$D_{p}\leq$') + str(psb_dub[psb_dubi]/1.e3) + ' '+ str(r'$\rm{\mu}$m')))
		

	ax0.set_ylabel(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	
	# include legend
	ax0.legend(fontsize=14)

	return()

def part_num_vs_time_sizeseg(self):

	import scipy.constants as si
	
	# number of actual particle size bins
	num_asb = (self.ro_obj.nsb-self.ro_obj.wf)

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))
	
	try: # in case user-supplied values given:
		psb_dub = [float(i) for i in self.e303p.text().split(',')]
	except:
		# default
		# particle size bins (upper bounds) diameter (um), 
		# ready for indexing the concentrations that will be summed
		psb_dub = [1.e-1, 5.e-1, 1.e0, 2.5e0, 1.e1, 2.e1]
	
	# convert diameters to nm from um
	for i in range(len(psb_dub)):
		psb_dub[i] = psb_dub[i]*1.e3

	# convert recorded radius to nm from um and from radius to diameter
	dbou_rec = self.ro_obj.rad*1.e3*2.
	
	# index array preparation
	psb_ind = np.zeros((len(self.ro_obj.thr), num_asb))

	# get the index of particle size bins
	for psb_dubi in psb_dub:
		psb_ind += dbou_rec[:, 1::]<psb_dubi

	# then use this index to find the mass concentrations 
	# inside each size bound

	# results array preparation
	psb_res = np.zeros((len(self.ro_obj.thr), len(psb_dub)))

	# temporary holder array
	nc_psb_temp = np.zeros((len(self.ro_obj.thr), num_asb))

	for psb_dubi in range(len(psb_dub)): # loop through size categories
		
		# temporary holder array
		nc_psb_temp[:, :] = self.ro_obj.Nrec_dry[:, :]

		# zero the unwanted size bins
		nc_psb_temp[psb_ind<(len(psb_dub)-psb_dubi)] = 0.

		# get the concentrations of components in this size category (# particles/cm3)
		psb_res[:, psb_dubi] = np.sum(nc_psb_temp, axis=1)
	
		ax0.plot(self.ro_obj.thr, psb_res[:, psb_dubi], label = str(str(r'$D_{p}$<') + str(psb_dub[psb_dubi]/1.e3) + ' '+ str(r'$\rm{\mu}$m')))
		

	ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
	ax0.set_ylabel(r'Concentration (# particles$\,$cm$\rm{^{-3}}$)', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	
	# include legend
	ax0.legend(fontsize=14)

	return()

def part_area_vs_time_sizeseg(self):

	import scipy.constants as si
		
	# get required variables from self
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	

	# number of actual particle size bins
	num_asb = (self.ro_obj.nsb-self.ro_obj.wf)

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))
	
	try: # in case user-supplied values given:
		psb_dub = [float(i) for i in self.e303p.text().split(',')]
	except:
		# default
		# particle size bins (upper bounds) diameter (um), 
		# ready for indexing the concentrations that will be summed
		psb_dub = [1.e-1, 5.e-1, 1.e0, 2.5e0, 1.e1, 2.e1]
	
	# convert diameters to nm from um
	for i in range(len(psb_dub)):
		psb_dub[i] = psb_dub[i]*1.e3

	# get particle-phase concentrations (# molecules/cm3)
	yp = yrec[:, self.ro_obj.nc:(self.ro_obj.nc*(num_asb+1))]
	
	# zero water
	yp[:, self.ro_obj.H2O_ind::self.ro_obj.nc] = 0.
	
	# convert # molecules/cm3 to moles/cm3
	yp = yp/si.N_A
	
	# get volume of every component (cm3/cm3)
	yp = yp*np.tile(np.array(self.ro_obj.comp_MV).reshape(1, -1), (1, num_asb))

	# repeat number concentration of particles over components (# particles/cm3)
	Nrec_wet_rep = np.repeat(self.ro_obj.Nrec_wet, self.ro_obj.nc, axis=1)

	# get volume of every component per particlecm3/cm3)
	yp[Nrec_wet_rep>0] = yp[Nrec_wet_rep>0]/Nrec_wet_rep[Nrec_wet_rep>0]
	yp[Nrec_wet_rep == 0] = 0.
	
	vp = np.zeros((yrec.shape[0], num_asb))
	
	# sum over components
	for i in range(num_asb):
		vp[:, i] = np.sum(yp[:, self.ro_obj.nc*i:self.ro_obj.nc*(i+1)], axis=1)
	
	# get radius of every size bin particle at every time step (cm/cm3)
	radius_dry = (((3.*vp)/(4.*np.pi))**(1./3.))

	# get surface area summed across all particles (cm2/cm3), assuming spherical shape
	sa = 4.*np.pi*radius_dry**2*self.ro_obj.Nrec_wet
	
	# convert surface area from cm2/cm3 to m2/m3
	sa = sa*1.e2

	# convert dry radius to nm from cm and from radius to diameter
	dbou_rec = radius_dry*1.e7*2.
	
	# index array preparation
	psb_ind = np.zeros((len(self.ro_obj.thr), num_asb))

	# get the index of particle size bins relative to user-defined size bins
	for psb_dubi in psb_dub:
		psb_ind += dbou_rec[:, :]<psb_dubi
		
	# then use this index to find the surface area concentrations 
	# inside each size bound

	# results array preparation
	psb_res = np.zeros((len(self.ro_obj.thr), len(psb_dub)))

	# temporary holder array
	sac_psb_temp = np.zeros((len(self.ro_obj.thr), num_asb))

	for psb_dubi in range(len(psb_dub)): # loop through size categories
		
		# temporary holder array
		sac_psb_temp[:, :] = sa

		# zero the unwanted size bins
		sac_psb_temp[psb_ind<(len(psb_dub)-psb_dubi)] = 0.

		# get the concentrations in this size category (m2/m3)
		psb_res[:, psb_dubi] = np.sum(sac_psb_temp, axis=1)
	
		ax0.plot(self.ro_obj.thr, psb_res[:, psb_dubi], label = str(str(r'$D_{p}$<') + str(psb_dub[psb_dubi]/1.e3) + ' '+ str(r'$\rm{\mu}$m')))
		
	ax0.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
	ax0.set_ylabel(r'Concentration (m$\rm{^{2}}\,$ m $\rm{^{-3}}$)', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	
	# include legend
	ax0.legend(fontsize=14)

	return()

# function for plotting contribution per user-supplied component by mass
# against time
def comp_part_mass_vs_time(self):

	import scipy.constants as si

	# inputs: ----
	# self - PyCHAM object
	# ------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], 
		self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], 
		self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], 
		self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = self.ro_obj.comp_MW
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], 
		self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]

	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))

	# components to consider
	cn = self.comp_names_part_mass_vs_time

	# prepare results matrix
	res = np.zeros((len(cn), len(timehr)))	

	# particle-phase concentrations of all 
	# components (# molecules/cm3)
	if (wall_on > 0): # wall on
		ppc = yrec[:, num_comp:-num_comp*wall_on]
	if (wall_on == 0): # wall off
		ppc = yrec[:, num_comp::]
	
	# tile molar masses over size bins and times
	y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, num_asb))
	y_mwt = np.tile(y_mwt, (len(timehr), 1))
	# convert from # molecules/cm3 to ug/m3
	ppc = (ppc/si.N_A)*y_mwt*1.e12

	# loop through components to plot
	for i in range(len(cn)):	

		# get chemical scheme index
		try:
			ci = comp_names.index(cn[i])
			# sum over size bins
			res[i, :] = np.sum(ppc[:, ci::num_comp], axis=1)
						
		except:
			if (cn[i] == 'SOA'):
				
				y_SOA = np.zeros((ppc.shape[0], 
					ppc.shape[1]))

				y_SOA[:, :] = ppc[:, :]

				# zero water and seed
				y_SOA[:, H2Oi::num_comp] = 0.
				# loop through seed components
				for si in seedi: 
					y_SOA[:, si::num_comp] = 0.
				
				# sum over components and size bins
				res[i, :] = np.sum(y_SOA, axis=1)
	
	ax0.stackplot(timehr[1:-1], res[:, 1:-1], labels = cn)
	plt.yscale('log')
	ax0.set_ylabel(r'Particle Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)', fontsize = 14)
	ax0.set_xlabel(r'Time through day (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14, loc='upper left')

	return()	

# function for plotting contribution per user-supplied component by
# risk (as a function of mass) against time
def comp_part_risk_vs_time(self):

	import scipy.constants as si

	# inputs: ----
	# self - PyCHAM object
	# ------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], 
		self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], 
		self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], 
		self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = self.ro_obj.comp_MW
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], 
		self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]

	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))

	# components to consider
	cn = self.comp_names_part_mass_vs_time

	# prepare results matrix
	res = np.zeros((len(cn), len(timehr)))	

	# particle-phase concentrations of all 
	# components (# molecules/cm3)
	if (wall_on > 0): # wall on
		ppc = yrec[:, num_comp:-num_comp*wall_on]
	if (wall_on == 0): # wall off
		ppc = yrec[:, num_comp::]
	
	# tile molar masses over size bins and times
	y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, num_asb))
	y_mwt = np.tile(y_mwt, (len(timehr), 1))
	# convert from # molecules/cm3 to ug/m3
	ppc = (ppc/si.N_A)*y_mwt*1.e12

	# concentration-response function from 
	# https://assets.publishing.service.gov.uk/
	# media/623302158fa8f504aa780865/
	# COMEAP_Statement_on_PM2.5_mortality_quantification.pdf
	risk_coeff = 1.08

	# loop through components to plot
	for i in range(len(cn)):	

		# get chemical scheme index
		try:
			ci = comp_names.index(cn[i])
			# sum over size bins
			res[i, :] = risk_coeff**((
			np.sum(ppc[:, ci::num_comp], axis=1))/10.)-1.
						
		except:
			if (cn[i] == 'SOA'):
				
				y_SOA = np.zeros((ppc.shape[0], 
					ppc.shape[1]))

				y_SOA[:, :] = ppc[:, :]

				# zero water and seed
				y_SOA[:, H2Oi::num_comp] = 0.
				# loop through seed components
				for si in seedi: 
					y_SOA[:, si::num_comp] = 0.
				
				# sum over components and size bins
				res[i, :] = risk_coeff**((np.sum(
				y_SOA, axis=1))/10.)-1.
	
	ax0.stackplot(timehr[1:-1], res[:, 1:-1], labels = cn)

	ax0.set_ylabel(r'factor change in all-cause mortality', 
	fontsize = 14)
	ax0.set_xlabel(r'Time through day (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14)

	return()


# function for plotting contribution per user-supplied component by
# carbon oxidation state against time
def comp_part_cos_vs_time(self):

	import scipy.constants as si

	# inputs: ----
	# self - PyCHAM object
	# ------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], 
		self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], 
		self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	Ndry = np.zeros((self.ro_obj.Nrec_dry.shape[0], 
		self.ro_obj.Nrec_dry.shape[1]))
	Ndry[:, :] = self.ro_obj.Nrec_dry[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = self.ro_obj.comp_MW
	HC = np.array((self.ro_obj.HyC))
	
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], 
		self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]

	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize = (14, 7))

	# components to consider
	cn = self.comp_names_part_mass_vs_time
	
	# remember to include total carbon oxidation state
	#cn.append('total')
	
	# prepare results matrix
	res = np.zeros((len(cn), len(timehr)))	
	ros = np.zeros((len(cn), len(timehr)))	

	# particle-phase concentrations of all 
	# components (# molecules/cm3)
	if (wall_on > 0): # wall on
		ppc = yrec[:, num_comp:-num_comp*wall_on]
	if (wall_on == 0): # wall off
		ppc = yrec[:, num_comp::]
	
	# tile molar masses over size bins and times
	y_mwt = np.tile(np.array((y_MW)).reshape(1, -1), (1, num_asb))
	y_mwt = np.tile(y_mwt, (len(timehr), 1))
	# convert from # molecules/cm3 to ug/m3
	ppc = (ppc/si.N_A)*y_mwt*1.e12

	# now transform ppc matrix so that times in rows,
	# components in columns and size bins in third dimension
	ppc = ppc.reshape(len(timehr), num_comp, num_asb, order='F')

	# now sum over size bins
	ppc = np.sum(ppc, axis=2)

	# carbon oxidation state per component 
	cos = np.zeros((len(rel_SMILES)))
	for ci in range(len(rel_SMILES)): # loop through components
		# get carbon number
		Cn = (rel_SMILES[ci].count('c')+
			rel_SMILES[ci].count('C'))
		if (Cn == 0):
			cos[ci] = 0.
		else:
			
			cos[ci] = (((HC[0, ci]*Cn)*-1+
			rel_SMILES[ci].count('o')*2+
			rel_SMILES[ci].count('O')*2.+
			rel_SMILES[ci].count('oo')*-2+
			rel_SMILES[ci].count('oO')*-2+
			rel_SMILES[ci].count('Oo')*-2+
			rel_SMILES[ci].count('OO')*-2+
			rel_SMILES[ci].count('O[O]')*-3+
			rel_SMILES[ci].count('ON(=O)=O')*-5)/Cn)
			
	# tile carbon oxidation states over times
	cos = np.tile(cos.reshape(1, -1), (len(timehr), 1))	
	# get the ROS
	ros0 = cos*1.1+1.51
	# can't have negative reactive oxidation spaces	
	ros0[ros0<0.] = 0. 

	# estimate mass fractions
	# zero mass of water
	ppc[:, H2Oi] = 0.
	#ppc[:, comp_names.index('core')] = 0.
	#ppc[:, comp_names.index('core_ind')] = 0.
	
	# mass fraction of each component at each time
	mf = (ppc/(np.sum(ppc, axis=1).reshape(-1, 1)))

	# multiply by mass concentration fraction
	cos_frac = cos*mf
	ros_frac = ros0*mf

	# loop through components of known 
	# oxidative potential
	op_res = np.zeros((6, len(timehr)))
	# ammonium sulphate (outdoor core)
	op_res[0, :] = 180.*mf[:, comp_names.index('AMM_SUL')]
	
	
	# outdoor primary (pri_org)
	op_res[1, :] = 9.*mf[:, comp_names.index('pri_org')]
	# outdoor secondary (sec_org), e.g. Figure 6 of Zhang et al. 2022:
	# https://doi.org/10.5194/acp-22-1793-2022
	op_res[2, :] = 60.*(
	mf[:, comp_names.index('sec_org1')]+
	mf[:, comp_names.index('sec_org0')]+
	mf[:, comp_names.index('sec_org-1')]+
	mf[:, comp_names.index('sec_org-2')])

	# indoor elemental carbon (bcin)
	op_res[3, :] = 2.*mf[:, comp_names.index('bcin')]
	# indoor primary organic
	op_res[4, :] = 9.*mf[:, comp_names.index('pri_orgin')]

	# SOA using the factor from e.g. Figure 6 of Zhang et al. 2022:
	# https://doi.org/10.5194/acp-22-1793-2022, but note that a factor
	# 67 was used here previously but I can't remember the
	# reference
	# zero mass fractions of seed components
	mf[:, seedi] = 0.
	mf[:, H2Oi] = 0.
	op_res[5, :] = 60.*np.sum(mf, axis=1)

	# multiply by total dry mass of PM
	op_res = op_res*(np.sum(ppc, axis=1).reshape(1, -1))

	# loop through components to plot
	#for i in range(len(cn)):	
		
		# get chemical scheme index
		#try:
			#ci = comp_names.index(cn[i])
			# sum over size bins
			#res[i, :] = cos_frac[:, ci]
			#ros[i, :] = ros_frac[:, ci]
		#except:	
			#if (cn[i] == 'SOA'):
				
				#y_SOA = np.zeros((cos_frac.shape[0], 
				#	cos_frac.shape[1]))
				#yr_SOA = np.zeros((ros_frac.shape[0], 
				#	ros_frac.shape[1]))	

				#y_SOA[:, :] = cos_frac[:, :]
				#yr_SOA[:, :] = ros_frac[:, :]
				# zero water and seed
				#y_SOA[:, H2Oi] = 0.
				#yr_SOA[:, H2Oi] = 0.
				
				# loop through seed components
				#for si in seedi: 
				#	y_SOA[:, si] = 0.
				#	yr_SOA[:, si] = 0.
					
				# sum over components
				#res[i, :] = np.sum(y_SOA, axis=1)
				#ros[i, :] = np.sum(yr_SOA, axis=1)
	
			#if (cn[i] == 'total'):
				#res[i, :] = np.sum(cos_frac, axis=1)
				#ros[i, :] = np.sum((ros[0:i, :]),
				#	 axis=0)
		#ax0.plot(timehr[1:-1]+16., res[i, 1:-1], label=cn[i])

	# now that we have ROS results for all components, we can
	# compare the ROS with and without each component, and 
	# therefore estimate its factor change in ROS	
	# loop through components contributing to ROS
	#for i in range(len(cn)-1):	
		#ros[i, :] = (ros[i, :]/ros[-1, :])
		#ax1.plot(timehr[1:-1]+16., ros[i, 1:-1], label=cn[i])
	print(cn)
	ax0.stackplot(timehr[1:-1], op_res[:, 1:-1]*1.e-3, labels=['AMM_SUL', 'pri_org', 'SOAout', 'bcin', 'pri_orgin', 'SOAin'])
	#ax0.set_ylabel(str('''Carbon oxidation state weighted by ''' +
	#'''mass'''), fontsize = 14)
	ax0.set_ylabel(str('''$\mathrm{OP_{DTT}}$ $\mathrm{(nmol\, min^{-1})}$ weighted by ''' +
	'''mass'''), fontsize = 14)
	#ax1.set_ylabel('''ROS''')
	ax0.set_xlabel(r'Time through day (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	#ax1.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	#ax1.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14)

	return()
