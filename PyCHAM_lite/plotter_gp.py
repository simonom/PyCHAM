########################################################################
#								       #
# Copyright (C) 2018-2024					       #
# Simon O'Meara : simon.omeara@manchester.ac.uk			       #
#								       #
# All Rights Reserved.                                                 #
# This file is part of PyCHAM                                          #
#                                                                      #
# PyCHAM is free software: you can redistribute it and/or modify it    #
# under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or    #
# (at  your option) any later version.                                 #
#                                                                      #
# PyCHAM is distributed in the hope that it will be useful, but        #
# WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  #
# General Public License for more details.                             #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with PyCHAM.  If not, see <http://www.gnu.org/licenses/>.      #
#                                                                      #
########################################################################
'''plots results for the gas-phase temporal profiles of specified 
components'''
# simulation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
# for customised colormap
from matplotlib.colors import LinearSegmentedColormap
# set colormap tick labels to standard notation
import matplotlib.ticker as ticker 
import os
import retr_out
import numpy as np
import scipy.constants as si
import openbabel.pybel as pybel

# plotting gas-phase concentrations
def plotter(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0 for ug/m3 or 1 for ppb, 3 for # molecules/cm3) 
	# or tests (2) are the calling module
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
	Cfac = self.ro_obj.cfac
	group_indx = self.ro_obj.gi 
	comp0 = self.ro_obj.init_comp
	
	y_MW = np.array(y_MW) # convert to numpy array from list
	Cfac = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list
	
	# number of actual particle size bins
	num_asb = (num_sb-wall_on)	

	if (caller == 0 or caller == 1 or caller == 3 or caller == 4 or caller == 5 or caller == 6):
		# show results to screen and turn on interactive mode
		plt.ion() 
	
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	if (comp_names_to_plot): # if component names specified
	
		# gas-phase concentration sub-plot --------------------	
		ip_fail = 0 # start by assuming all requested components available

		for i in range(len(comp_names_to_plot)):

			group_flag = 0 # start by assuming single components wanted
			
			if (comp_names_to_plot[i].strip() == 'H2O'):
				indx_plot = [H2Oi]
				indx_plot = np.array((indx_plot))
			if (comp_names_to_plot[i].strip() == 'RO2'):
				indx_plot = (np.array((group_indx['RO2i'])))
				group_flag = 1
			if (comp_names_to_plot[i].strip() == 'RO'):
				indx_plot = (np.array((group_indx['ROi'])))
				group_flag = 1	
			if (comp_names_to_plot[i].strip() == 'HOMRO2'):
				indx_plot = (np.array((group_indx['HOMRO2'])))
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1	
			if (comp_names_to_plot[i].strip() == 'HOM'):
				indx_plot = (np.array((group_indx['HOMs'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1	
			if (comp_names_to_plot[i].strip() == '-OOH'):
				indx_plot = (np.array((group_indx['OOH'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			if (comp_names_to_plot[i].strip() == 'HOM-OOH'):
				indx_plot = (np.array((group_indx['HOM_OOH'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1	
			if (comp_names_to_plot[i].strip() == '-OH'):
				indx_plot = (np.array((group_indx['OH'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1			
			if (comp_names_to_plot[i].strip() == 'HOM-OH'):
				indx_plot = (np.array((group_indx['HOM_OH'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			if (comp_names_to_plot[i].strip() == '-carbonyl'):	
				indx_plot = (np.array((group_indx['carbonyl'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			if (comp_names_to_plot[i].strip() == 'HOM-carbonyl'):
				indx_plot = (np.array((group_indx['HOM_carbonyl'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			if (comp_names_to_plot[i].strip() == '-NO3'):
				indx_plot = (np.array((group_indx['NO3'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1
			if (comp_names_to_plot[i].strip() == 'HOM-NO3'):
				indx_plot = (np.array((group_indx['HOM_NO3'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1	
			if (comp_names_to_plot[i].strip() == 'ROOR'):
				indx_plot = (np.array((group_indx['ROOR'])))			
				group_flag = 1
				if (indx_plot.shape[0] == 0):
					ip_fail = 1	

			if ('C' in comp_names_to_plot[i].strip() or 'c' in comp_names_to_plot[i].strip()):
				# if a number given after carbon atom
				if comp_names_to_plot[i].strip()[1::].isnumeric():
					# get carbon number
					Cn = float(comp_names_to_plot[i].strip()[1::])
					# get all components with this many carbons
					# empty list
					indx_plot = []
					indx_cnt = 0 # keep count on species
					for SMILESi in rel_SMILES:
						if (SMILESi.count('C') + SMILESi.count('c')) == Cn:
							indx_plot.append(indx_cnt)
						indx_cnt += 1 # keep count on species

					indx_plot = (np.array((indx_plot)))			
					group_flag = 1
					if (indx_plot.shape[0] == 0):
						ip_fail = 1
			if (ip_fail == 1):
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

			if (comp_names_to_plot[i].strip() != 'H2O' and group_flag != 1):
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
			
			if (caller == 0  or caller == 5): # ug/m3 plot
			
				# gas-phase concentration (# molecules/cm3)
				conc = yrec[:, indx_plot].reshape(yrec.shape[0], (indx_plot).shape[0])*Cfac

				# gas-phase concentration (ug/m3)
				conc = ((conc/si.N_A)*y_MW[indx_plot])*1.e12
			
			if (caller == 1 or caller == 4): # ppb plot
			
				# gas-phase concentration (ppb)
				conc = yrec[:, indx_plot].reshape(yrec.shape[0], (indx_plot).shape[0])

			if (caller == 3  or caller == 6): # # molecules/cm3 plot
			
				# gas-phase concentration (# molecules/cm3)
				conc = yrec[:, indx_plot].reshape(yrec.shape[0], (indx_plot).shape[0])*Cfac
			
			if (len(indx_plot) > 1): # e.g. for sum of RO2 or RO
				conc = np.sum(conc, axis=1) # sum multiple components
			
			# plot this component
			if (group_flag == 0): # if not a sum of components
				if (caller == 4 or caller == 5 or caller == 6): # log10 y axis
					ax0.semilogy(timehr, conc, '-+', linewidth = 4., label = str(str(comp_names[int(indx_plot)]+' (gas-phase)')))
				if (caller == 0 or caller == 1 or caller == 3): # linear y axis
					ax0.plot(timehr, conc, '-+', linewidth = 4., label = str(str(comp_names[int(indx_plot)]+' (gas-phase)')))
			
			else: # if a sum over a group of components
				if (caller == 4 or caller == 5 or 
					caller == 6): # log y axis
					ax0.semilogy(timehr, conc, '-+',
					 linewidth = 4., label = 
					str(r'$\Sigma$' + 
					comp_names_to_plot[i].strip() + 
					' (gas-phase)'))
				if (caller == 0 or caller == 1 or caller == 3): # linear y axis
					ax0.plot(timehr, conc, '-+', linewidth = 4., label = str(r'$\Sigma$' + comp_names_to_plot[i].strip() + ' (gas-phase)'))

		if (caller == 0 or caller == 5): # ug/m3 plot
			ax0.set_ylabel(str(r'Concentration ($\rm{\mu}$g$\,$m$\rm{^{-3}}$)'), fontsize = 14)
		if (caller == 1 or caller == 4): # ppb plot
			ax0.set_ylabel(str(r'Mixing ratio (ppb)'), fontsize = 14)
		if (caller == 3 or caller == 6): # # molecules/cm3 plot
			gpunit = str('\n(' + u'\u0023' + ' molecules/cm' + u'\u00B3' + ')')
			ax0.set_ylabel(str(r'Concentration ' + gpunit), fontsize = 14)

		ax0.set_xlabel(str(r'Time through simulation (hours)'), fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.legend(fontsize = 14)
		
		# end of gas-phase concentration sub-plot --------------
	

	if (caller == 2): # display
		plt.show()	

	return()


# for time series of the average organic peroxy radical molecule
def RO2_av_molec(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------
	# caller - marker for whether PyCHAM (0 for ug/m3 or 1 for ppb, 3 for # molecules/cm3) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------

	# chamber condition ---------------------------------------------------------
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	
	y_MW = np.array(y_MW) # convert to numpy array from list
	Cfac = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	# get RO2 indices		
	indx_plot = (np.array((group_indx['RO2i'])))

	# gas-phase concentration (# molecules/cm3)
	conc = yrec[:, indx_plot].reshape(yrec.shape[0], (indx_plot).shape[0])*Cfac

	# sum of concentrations (# molecules/cm3)
	conc_sum = np.sum(conc, axis = 1)

	# get SMILES of RO2 componenets
	rel_SMILES = rel_SMILES[indx_plot]

	# number of carbons and oxygens in each component
	Ccnt = []
	Ocnt = []
	Hcnt = []
	for i in rel_SMILES:
		Ccnt.append(i.count('C')+i.count('c'))
		Ocnt.append(i.count('O'))
		# generate pybel object
		Pybel_object = pybel.readstring('smi', i)
		try:
			Hi = (Pybel_objects[indx].formula).index('H')
			Hcnt = -1
		except:
			Hcnt = 0.0
		if (Hcnt == -1):
			try:
				Hcnt.append(float(Pybel_objects[indx].formula[Hi+1:Hi+3]))
			except:
				Hcnt.append(float(Pybel_objects[indx].formula[Hi+1:Hi+2]))
	Ccnt = np.tile(((np.array((Ccnt))).reshape(1, -1)), (conc.shape[0], 1))
	Ocnt = np.tile(((np.array((Ocnt))).reshape(1, -1)), (conc.shape[0], 1))
	Hcnt = np.tile(((np.array((Hcnt))).reshape(1, -1)), (conc.shape[0], 1))
	
	# average carbon number of organic peroxy radicals (RO2) in gas-phase at each time step
	Cav_cnt = ((np.sum(Ccnt*conc, axis = 1))/conc_sum)
	# average oxygen number of organic peroxy radicals (RO2) in gas-phase at each time step
	Oav_cnt = ((np.sum(Ocnt*conc, axis = 1))/conc_sum)
	# average hydrogen number of organic peroxy radicals (RO2) in gas-phase at each time step
	Hav_cnt = ((np.sum(Hcnt*conc, axis = 1))/conc_sum)

	ax0.plot(timehr, Cav_cnt, '-+', linewidth = 4., label = 'Carbon number')
	ax0.plot(timehr, Oav_cnt, '-+', linewidth = 4., label = 'Oxygen number')
	ax0.plot(timehr, Hav_cnt, '-+', linewidth = 4., label = 'Hydrogen number')
	
	ax0.set_ylabel(str(r'Average number of atoms per RO2 molecule'), fontsize = 14)
	ax0.set_xlabel(str(r'Time through simulation (hours)'), fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14)

def plotter_noncsv(caller, dir_path, comp_names_to_plot, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0 for ug/m3 or 1 for ppb, 3 for # molecules/cm3) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# comp_names_to_plot - chemical scheme names of components to plot
	# self - reference to GUI
	# --------------------------------------------------------------------------
	
	# get required information from EASY
	[Etime_s, Ecomp_names, ECrec, TEMP, PRESS] = retr_out.retr_out_noncsv(dir_path, comp_of_int)
	
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	TEMP = 293.15 # chamber temperature (K)
	PRESS = 1.e5 # chamber pressure (Pa)
	ntot = PRESS*(si.N_A/((si.R*1.e6)*TEMP))
	# one billionth of number of molecules in chamber unit volume
	Cfac = (ntot*1.e-9) # ppb to molecules/cm3 conversion factor
	
	if (comp_names_to_plot): # if component names specified
	
		# gas-phase concentration sub-plot ---------------------------------------------	
		for i in range(len(comp_names_to_plot)):
			
			Ei = Ecomp_names.index(comp_names_to_plot[i]) # EASY index
			ax0.plot(Etime_s/3600., (ECrec[:, Ei]/Cfac), '-x', linewidth = 2., label = str(comp_names_to_plot[i]))
		
	
	ax0.set_ylabel(str(r'Concentration (ppb)'), fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.legend(fontsize = 14)
	if (caller <= 2): # display
		plt.show()	
	
	return()

# plotting the radical pool
def plotter_rad_pool(self):

	# ---------------------
	# inputs:
	# self - reference to PyCHAM
	# ---------------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	group_indx = self.ro_obj.gi
	OC = self.ro_obj.O_to_C
	PsatPa = self.ro_obj.vpPa
	gen_num = np.zeros((len(self.ro_obj.gen_numbers), 1))
	gen_num[:, 0] = self.ro_obj.gen_numbers
	
	if (self.rad_mark == 0): # if alkyl peroxy radicals
		# get RO2 indices
		indx_plot = np.array((group_indx['RO2i']))
		
	if (self.rad_mark == 1): # if alkoxy radicals
		# get RO indices		
		indx_plot = np.array((group_indx['ROi']))
	
 	# need a section here to conditionally 
	# ignore radicals with carbon number below the user-supplied
	comp_cnt = 0
	for SMILEi in rel_SMILES: # loop through SMILE strings
		# if this component has less than the threshold carbon number
		if (SMILEi.count('C')+SMILEi.count('c') < self.Cnum_thresh):
			if comp_cnt in indx_plot:
				indx_plot = np.delete(indx_plot, indx_plot == comp_cnt)
		comp_cnt += 1		
	
	# get names of radicals in this pool
	rad_names = (np.array((comp_names)))[indx_plot]

	# isolate gas-phase concentrations of radical (ppb)
	y_rad = yrec[:, indx_plot].reshape(yrec.shape[0], (indx_plot).shape[0])
	
	# isolate generations of radicals
	gen_num = gen_num[indx_plot, 0]

	# prepare for fractional contribution
	y_radf = np.zeros((y_rad.shape[0], y_rad.shape[1]))

	# sum of contributions per time step
	rad_sum = ((np.sum(y_rad, axis= 1)).reshape(-1, 1))

	# fractional contribution per time step
	y_radf = np.zeros((y_rad.shape[0], y_rad.shape[1]))
	y_radf[rad_sum[:,0]>0, :] = y_rad[rad_sum[:,0]>0, :]/rad_sum[rad_sum[:,0]>0, :]
	
	# sum of fractional contributions over time per component
	y_radf_tot = np.sum(y_radf, axis = 0)

	# order of radicals with greatest contributor first
	ord = y_radf_tot.argsort()
	
	# use just the indices of the top number given by user
	ord = ord[-self.rad_ord_num::]
	
	# plot fractional contribution ---------------------
	plt.ion() # show results to screen and turn on interactive mode
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	par1 = ax0.twinx() # parasite (right) axis	

	frac_sum = np.zeros((len(timehr))) # results array for total fractions shown here
	
	for i in range(len(ord)): # loop through top contributors, with biggest first
		ax0.plot(timehr, y_radf[:, ord[-(i+1)]], label = str(rad_names[ord[-(i+1)]] + ' (gen' + str(int(gen_num[ord[-(i+1)]])) + ')' + ' frac.'))
		frac_sum += y_radf[:, ord[-(i+1)]] # total fractions shown here
		# plot right axis (absolute concentration)
		p3, = par1.plot(timehr, y_rad[:, ord[-(i+1)]], '--', label = str(rad_names[ord[-(i+1)]] + ' conc.'))

	# also plot sum of fractions shown in plot
	ax0.plot(timehr, frac_sum, '-k', 
		label = str(r'$\Sigma$(frac. shown here)'))

	# in case you want to check that sum of fractions=1
	#ax0.plot(timehr, np.sum(y_radf, axis=1), 
	#	label = 'sum of fractions (check)')

	ax0.set_ylabel(str('Fraction of all concentrations'), 
		fontsize = 14)
	# right vertical axis label
	par1.set_ylabel(str(r'Concentration ($\mathrm{ppb}$)'), 
		fontsize = 14, rotation=270, labelpad=20) 
	ax0.set_xlabel(str(r'Time through simulation (hours)'), 
		fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	par1.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e')) # set tick format for vertical axis
		
	ax0.legend(fontsize = 14, loc= 'upper left')
	par1.legend(fontsize = 14, loc= 'upper right')

	return()

# plotting the radical flux
def plotter_rad_flux(self):

	# ---------------------
	# inputs:
	# self - reference to PyCHAM
	# ---------------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	
	# no record of change tendency for final experiment time point
	timehr = timehr[0:-1]

	if (self.rad_mark == 2): # if alkyl peroxy radicals
		# get RO2 indices		
		indx_plot = np.array((group_indx['RO2i']))

	if (self.rad_mark == 3): # if alkoxy radicals
		# get RO indices		
		indx_plot = np.array((group_indx['ROi']))
	
	# get names of radicals in this pool
	rad_names = (np.array((comp_names)))[indx_plot]

	# check that all these components present
	for comp_name in (rad_names):
		
		fname = str(self.dir_path+ '/' + comp_name + '_rate_of_change')
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
	par1 = ax0.twinx() # parasite (right) axis	


	# empty results for sum of absolute values for production and loss terms
	cr_dydt = np.zeros((len(timehr), len(rad_names)))
	
	compi = 0 # count on components
	for comp_name in (rad_names): # loop through possible components to plot
		
		# open results for this component
		fname = str(self.dir_path+ '/' + comp_name + '_rate_of_change')
		dydt = np.loadtxt(fname, delimiter = ',', skiprows = 1) # skiprows = 1 omits header	
		
		ci = comp_names.index(comp_name) # get index of this component
		
		# note that penultimate column in dydt is gas-particle 
		# partitioning and final column is gas-wall partitioning, whilst
		# the first row contains chemical reaction numbers
		
		# sum chemical reaction gains
		crg = np.zeros((dydt.shape[0]-1, 1))
		# sum chemical reaction losses
		crl = np.zeros((dydt.shape[0]-1, 1))

		for ti in range(dydt.shape[0]-1): # loop through times
			indx = dydt[ti+1, 0:-2] > 0 # indices of reactions that produce component
			cr_dydt[ti, compi] = dydt[ti+1, 0:-2][indx].sum()
			indx = dydt[ti+1, 0:-2] < 0 # indices of reactions that lose component
			cr_dydt[ti, compi] += np.abs(dydt[ti+1, 0:-2][indx].sum())
		
		compi += 1 # count on components
	
	# sum fluxes over components
	rad_sum = (np.sum(cr_dydt, axis = 1)).reshape(-1, 1)

	# fraction of flux per time step
	cr_dydt_frac = np.zeros((cr_dydt.shape[0], cr_dydt.shape[1]))
	cr_dydt_frac[rad_sum[:, 0] > 0, :] = cr_dydt[rad_sum[:, 0] > 0, :]/rad_sum[rad_sum[:, 0] > 0, :]

	# sum fluxes per component across all times
	sum_cr_dydt = np.sum(cr_dydt, axis=0)
	# get order
	ord = sum_cr_dydt.argsort()
	
	# loop through top contributors
	for compi in range(self.rad_ord_num):
		
		# plot temporal profiles of fractional change 
		# tendencies due to chemical 
		# reaction
		ax0.plot(timehr, cr_dydt_frac[:, ord[-(compi+1)]], 
			label = str(rad_names[ord[-(compi+1)]] + 
				' frac.'))
		# plot right axis (absolute flux)
		p3, = par1.plot(timehr, cr_dydt[:, ord[-(compi+1)]],
			 '--', label = str(rad_names[ord[-(compi+1)]] +
			 ' flux'))


	ax0.set_xlabel('Time through experiment (hours)', fontsize = 14)
	ax0.set_ylabel('Fraction of all change tendencies', 
		fontsize = 14)
	par1.set_ylabel(str(r'''Change tendency ($\mathrm{molecules\,
		cm^{-3}\, s^{-1}}$)'''), fontsize = 14, rotation=270,
		 labelpad=20) # right vertical axis label
		

	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	par1.yaxis.set_tick_params(labelsize = 14, direction = 'in')

	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
		
	ax0.legend(fontsize = 14, loc= 'upper left')
	par1.legend(fontsize = 14, loc= 'upper right')
			 

	return()

# function to plots ozone isopleths -----------------------------------------
def O3_iso(self):

	import pickle # for dealing with pickle files

	# ---------------------
	# inputs:
	# self - reference to PyCHAM
	# ---------------------

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = np.zeros((self.ro_obj.Nrec_wet.shape[0], self.ro_obj.Nrec_wet.shape[1]))
	Nwet[:, :] = self.ro_obj.Nrec_wet[:, :]
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	gen_num = np.zeros((len(self.ro_obj.gen_numbers), 1))
	gen_num[:, 0] = self.ro_obj.gen_numbers
	Cfac = self.ro_obj.cfac
	Cfac = np.squeeze((np.array(Cfac))) # convert to numpy array from list

	# get SMILE strings of all zero generation components,
	# ignoring seed and water
	zg_smile = np.array(rel_SMILES[0:-2])[gen_num[:, 0] == 0]
	# sum carbon numbers in SMILE strings
	for zg_i in zg_smile:
		if (zg_i.count('c') + zg_i.count('C') > 1): # more than 1 carbon
			if (zg_i.count('o') + zg_i.count('O') == 0): # no oxygen
				VOC_smile = zg_i # get VOC SMILE string
				break

	# VOC index
	self.VOCi = rel_SMILES.index(VOC_smile)
	
	# get range (# molecules/cm3) of gas-phase VOC seen in simulation
	VOC_range = [min(yrec[:, self.VOCi]*Cfac), max(yrec[:, self.VOCi]*Cfac)]

	# NOx index
	self.NOxi = [comp_names.index('NO'), comp_names.index('NO2'), comp_names.index('NO3')]
	self.NOi = comp_names.index('NO')
	self.NO2i = comp_names.index('NO2')

	# O3 index
	self.O3i = comp_names.index('O3')

	NOxsum = np.zeros((len(timehr)))
	# sum NOx constituents (# molecules/cm3)
	for NOxii in self.NOxi:
		NOxsum += yrec[:, NOxii]*Cfac
	# range of NOx (# molecules/cm3)
	NOx_range = [min(NOxsum), max(NOxsum)]
	NOx_range = [1.0*Cfac[0], 1500.*Cfac[0]]
	VOC_range = [1.0*Cfac[0], 200.*Cfac[0]]
	
	NOx_values = np.arange(NOx_range[0], NOx_range[1]*1.01, (NOx_range[1]-NOx_range[0])/3.)
	VOC_values = np.arange(VOC_range[0], VOC_range[1]*1.01, (VOC_range[1]-VOC_range[0])/3.)
	
	# using the initial conditions specified in the 
	# model variables file, but with NOx and VOC 
	# fixed at the values set here, identify the
	# equilibrium O3 concentration
	
	# empty results (VOC changes by row, NOx by column)
	O3_res = np.zeros((len(VOC_values), len(NOx_values)))

	Nc = 0 # NOx count
	for NOxvi in NOx_values: # loop through NOx values
		Vc = 0 # VOC count
		for VOCvi in VOC_values: # loop through VOC values
			
			# load original self attributes and values
			with open(os.path.join(self.dir_path, 'simulation_self.pickle'), "rb") as f:
				self_dict = pickle.load(f)
			
			# transfer original values to self for this run
			for key in self_dict:# loop through dictionary entries here
				# set object's attribute value here
				# if original entry tagged as original, then transform to functional
				if key[-5::] == '_orig':
					setattr(self, key[0:-5], self_dict[key])
				else:
					setattr(self, key, self_dict[key])
			
			# get equilibrium O3 concentration
			self.testf = 5 # modify test flag value
		
			# set NOx and VOC value
			self.NOxequil = NOxvi
			self.VOCequil = VOCvi

			# ensure time step sufficiently long for spin-up
			self.update_stp = 6.e1
			self.save_step = 1.8e4
			self.tot_time = 1.8e4
			# now run program
			from middle import middle # prepare to communicate with main program
		
			note_messf = 0 # cancel note message flag
		
			for prog in middle(self): # call on modules to solve problem
		
				if (isinstance(prog, str)): # check if it's a message
					mess = prog
					if (mess[0:4] == 'Stop'): # if it's an error message
						return()
			# O3 result (ppb)
			O3_res[Vc, Nc] = self.O3equil/Cfac[0]
			
			Vc += 1 # VOC count
		Nc += 1 # NOx count
	
	# prepare plot
	plt.ion() # display figure in interactive mode
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	# ozone contour plot (# molecules/cm3)
	p1 = plt.contourf(NOx_values/Cfac[0], VOC_values/Cfac[0], O3_res)
	# format axis labels
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')

	# set axis titles
	ax0.set_xlabel(str(r'NOx mixing ratio (ppb)'), fontsize=14)
	ax0.set_ylabel(str(r'VOC mixing ratio (ppb)'), fontsize=14)

	# colour bar
	cb = plt.colorbar(p1, pad=0.25, ax=ax0)
	# colour bar label
	cb.set_label('O3 mixing ratio (ppb)', size=14, rotation=270, labelpad=20)
	# format color bar label
	cb.ax.tick_params(labelsize=14)

	return()
