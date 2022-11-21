##########################################################################################
#                                                                                        											 #
#    Copyright (C) 2018-2022 Simon O'Meara : simon.omeara@manchester.ac.uk                  				 #
#                                                                                       											 #
#    All Rights Reserved.                                                                									 #
#    This file is part of PyCHAM                                                         									 #
#                                                                                        											 #
#    PyCHAM is free software: you can redistribute it and/or modify it under              						 #
#    the terms of the GNU General Public License as published by the Free Software       					 #
#    Foundation, either version 3 of the License, or (at your option) any later          						 #
#    version.                                                                            										 #
#                                                                                        											 #
#    PyCHAM is distributed in the hope that it will be useful, but WITHOUT                						 #
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS       			 #
#    FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more              				 #
#    details.                                                                            										 #
#                                                                                        											 #
#    You should have received a copy of the GNU General Public License along with        					 #
#    PyCHAM.  If not, see <http://www.gnu.org/licenses/>.                                 							 #
#                                                                                        											 #
##########################################################################################
'''plots results from model and observations'''
# simulation and observation results are represented graphically

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import matplotlib.cm as cm
import os
import numpy as np
import scipy.constants as si
import pybel
import openpyxl # for opening excel file

def plotter_gp_mod_n_obs(self): # for gas-phase concentration temporal profiles
	
	# open xls file
	wb = openpyxl.load_workbook(filename = self.xls_path)
	obs = wb['PyCHAM'] # get first sheet
	# empty list to hold information on handling observation data
	obs_setup = []
	# loop through columns to get information required for 
	# handling observational data
	for co in obs.iter_cols(values_only=True):
		obs_setup.append(co[0])
	
	# open worksheet with observations
	obs = wb[obs_setup[0]]
	xrs = int(obs_setup[1].split(':')[0])-1 # starting row for x-axis
	xre = int(obs_setup[1].split(':')[1])-1 # finishing row for x-axis
	xcol = int(obs_setup[2].split(':')[1])-1 # x-axis column
	yrs = int(obs_setup[6].split(':')[0])-1 # starting row for y-axis
	yre = int(obs_setup[6].split(':')[1])-1 # finishing row for y-axis
	ycs = int(obs_setup[7].split(':')[0])-1 # starting column for y-axis
	yce = int(obs_setup[7].split(':')[1])-1 # finishing column for y-axis
	obsx = np.zeros((xre-xrs)) # empty array for x-axis data
	obsy = np.zeros((yre-yrs, yce-ycs+1)) # empty array for y-axis data
	cn = 0 # index counter
	xcn = 0 # x-axis counter
	ycn = 0 # y-axis counter
	# get x-axis data and y-axis data
	for co in obs.iter_rows(values_only=True):
		if cn > xrs:
			obsx[xcn] = co[xcol]
			obsy[ycn, :] = co[ycs:yce+1]
			xcn += 1 # x-axis counter
			ycn += 1 # y-axis counter
		cn += 1 # index counter
		if cn> xre:
			break # stop loop through rows
	
	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	
	# legend labels of observed components
	llab = [lli for lli in obs_setup[9].split(',')]

	# loop through observed components to plot in order to plot
	for i in range(yce-ycs+1):
		if (self.gp_units[-4::] == 'near'): # linear y-axis
			eby = obsy[:, i]*0.20 # error bar array
			markers, caps, bars = ax0.errorbar(obsx, obsy[:, i], yerr = eby, label = llab[i])
			# loop through error bars to set transparency
			[bar.set_alpha(0.1) for bar in bars]
		if (self.gp_units[-4::] == 'log.'): # logarithmic y-axis
			ax0.semilogy(obsx, obsy[:, i], label = llab[i])
	
	# ----------------------------------------------------------
	# deal with model results
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = self.ro_obj.yrec
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = self.ro_obj.Nrec_wet
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = self.ro_obj.rad
	space_mode = self.ro_obj.spacing

	# subtract any time before lights on
	timehr += float(obs_setup[4])

	# loop through modelled components
	for mci in range(len(self.gp_names)):
		if (self.gp_names[mci].strip() == ''): # if empty
			continue
		try:
			indx_plt = comp_names.index(self.gp_names[mci].strip())
		except:
			self.l203a.setText(str('Component ' + self.gp_names[mci].strip() + ' not found in chemical scheme used for this simulation'))
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
			
		if (self.gp_units[-4::] == 'near'): # linear y-axis
			ax0.plot(timehr, yrec[:, indx_plt], '--', label = str(self.gp_names[mci] + ' sim.'))
		if (self.gp_units[-4::] == 'log.'): # logarithmic y-axis
			ax0.semilogy(timehr, yrec[:, indx_plt], '--', label = str(self.gp_names[mci] + ' sim.'))
	
	# x-axis title
	ax0.set_xlabel(obs_setup[3], fontsize = 14)
	# y-axis title
	ax0.set_ylabel(obs_setup[8], fontsize = 14)
	# set label font size
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')
	# include legend
	ax0.legend(fontsize = 14)

	# -----------------------------------------------------------

	return()

def plotter_VK_mod_n_obs(self): # for Van Krevelen diagrams
	
	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = self.ro_obj.yrec
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = self.ro_obj.Nrec_wet
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = self.ro_obj.rad
	space_mode = self.ro_obj.spacing
	
	self.HC = np.array((ro_obj.HyC)) # hydrogen:carbon ratios of each component


	if (self.oandm == 2 or self.oandm == 4 or self.oandm == 6): # if user wants average over Hydrocarbons
		# index of the required experiment period
		tindx = (timehr*3600. >= 0)*(timehr*3600. <= timehr[-1]*3600.)
		# isolate gas-phase concentrations during this period (ppb)
		yrel = yrec[tindx, 0:num_comp]
		
		if (self.oandm == 2): # all hydrocarbons wanted
			# index for all hydrocarbons
			HCindx = (self.HC > 0.)
		
		# non-methane hydrocarbons wanted
		if (self.oandm == 4): # non-methane hydrocarbons wanted
			
			Ccount = [] # count number of carbons in each component
			for comp in rel_SMILES:
				Ccount.append(comp.count('C')+comp.count('c'))

			# index for non-methane hydrocarbons
			HCindx = (self.HC > 0.)*(np.array((Ccount)) > 1)

		# _ extension hydrocarbons wanted
		if (self.oandm == 6): # _ extension hydrocarbons wanted
			
			ext_cnt = [] # count number of carbons in each component
			for comp in comp_names:
				if comp.count('_') > 0:
					ext_cnt.append(1)
				else:
					ext_cnt.append(0)

			# index for _ extension hydrocarbons
			HCindx = (np.array((ext_cnt)) > 0)

		# isolate hydrocarbon molecules
		yrel = yrel[:, HCindx]
		# sum of concentrations
		HCsum = yrel.sum(axis=1).reshape(-1, 1)

		# fractional gas-phase concentration
		z_indx = (HCsum > 0.).reshape(-1, 1)

		yrel_frac = yrel[z_indx[:, 0], :]/((HCsum[z_indx]).reshape(-1, 1))
		OC_rel = (np.array((OC)))[HCindx].reshape(1, -1) # and ensure correct alignment with yrel
		HC_rel = self.HC[HCindx].reshape(1, -1) # and ensure correct alignment with yrel
		# weight concentrations by O:C ratio
		OCweight = yrel_frac*OC_rel
		# weight concentrations by H:C ratio
		HCweight = yrel_frac*HC_rel
	
		# get the average gas-phase O:C and H:C ratios during the stated period
		OCav = OCweight.sum(axis=1)
		HCav = HCweight.sum(axis=1)

		plt.ion() # show results to screen and turn on interactive mode
		
		# prepare plot
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

		# colormap from library (https://matplotlib.org/stable/api/cm_api.html#matplotlib.cm.get_cmap)
		colors = cm.get_cmap('ocean')

		# set contour levels
		levels = (MaxNLocator(nbins = 256).tick_values(np.min(timehr[tindx]*3600), np.max(timehr[tindx]*3600)))
	
		# associate colours and contour levels
		norm1 = BoundaryNorm(levels, ncolors=256, clip=True)
	
		for i in range(len(OCav)):
			ax0.plot(OCav[i], HCav[i], 'o', mec = 'k', mfc = colors(i/(len(OCav)-1)))

		# main title
		if (self.oandm == 2):
			ax0.set_title('All Hydrocarbons')
		if (self.oandm == 4):
			ax0.set_title('Non-methane Hydrocarbons')
		if (self.oandm == 6):
			ax0.set_title('Hydrocarbons Containing _ in Name')
		# x-axis title
		ax0.set_xlabel('O:C ratio number-averaged over hydrocarbons', fontsize = 14)
		# y-axis title
		ax0.set_ylabel('H:C ratio number-averaged over hydrocarbons', fontsize = 14)
		# set label font size
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

		# include colorbar to demonstrate time through experiment (s)
		cb = fig.colorbar(cm.ScalarMappable(norm=norm1, cmap=colors), ax=ax0)

		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('Time since experiment start (s)', size=14, rotation=270, labelpad=20)

	if (self.oandm == 3 or self.oandm == 5 or self.oandm == 7): # if user wants individual hydrocarbons at a given time

		# get the time through experiment (s) in question
		try:
			treq = float((self.e403.toPlainText()))
		except:
			treq = timehr[-1]*3600.

		# index of the required experiment period
		tindx = np.min(np.abs(timehr*3600-treq)) == (np.abs(timehr*3600-treq))

		# isolate gas-phase concentrations during this period (ppb)
		yrel = yrec[tindx, 0:num_comp]

		if (self.oandm == 3): # all individual hydrocarbons
			HCindx = self.HC > 0.

		if (self.oandm == 5): # non-methane individual hydrocarbons
			Ccount = [] # count number of carbons in each component
			for comp in rel_SMILES:
				Ccount.append(comp.count('C')+comp.count('c'))

			# index for non-methane hydrocarbons
			HCindx = (self.HC > 0.)*(np.array((Ccount)) > 1)

				# _ extension hydrocarbons wanted

		if (self.oandm == 7): # _ extension hydrocarbons wanted
			
			ext_cnt = [] # count number of carbons in each component
			for comp in comp_names:
				if comp.count('_') > 0:
					ext_cnt.append(1)
				else:
					ext_cnt.append(0)

			# index for _ extension hydrocarbons
			HCindx = (np.array((ext_cnt)) > 0)

		# isolate hydrocarbon molecules
		yrel = yrel[0, HCindx]

		# isolate H:C and O:C of hydrocarbons
		OC = (np.array(OC))[HCindx]
		HC = (np.array(self.HC))[HCindx]
		comp_names = (np.array(comp_names))[HCindx]
		
		# sum of these hydrocarbons
		HC_sum = np.sum(yrel)

		if (HC_sum > 0.):
			# normalise (0-1) hydrocarbon mixing ratios
			yrel = yrel/HC_sum
		else:
			yrel = []
			OC = []
			HC = []

		plt.ion() # show results to screen and turn on interactive mode
		
		# prepare plot
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

		# colormap from library (https://matplotlib.org/stable/api/cm_api.html#matplotlib.cm.get_cmap)
		colors = cm.get_cmap('ocean')

		# set contour levels
		levels = (MaxNLocator(nbins = 256).tick_values(0., 1.))
	
		# associate colours and contour levels
		norm1 = BoundaryNorm(levels, ncolors=256, clip=True)
	
		for i in range(len(HC)):
			ax0.plot(OC[i], HC[i], 'o', mec = 'k', mfc = colors(yrel[i]))
			
		# plot title
		if (self.oandm == 3):
			ax0.set_title(str('Van Krevelen at ' + str(treq) + ' s since experiment start for all hydrocarbons'), fontsize = 14)
		if (self.oandm == 5):
			ax0.set_title(str('Van Krevelen at ' + str(treq) + ' s since experiment start for non-methane hydrocarbons'), fontsize = 14)
		if (self.oandm == 7):
			ax0.set_title(str('Van Krevelen at ' + str(treq) + ' s since experiment start for extension hydrocarbons'), fontsize = 14)
		
		# x-axis title
		ax0.set_xlabel('O:C ratio of hydrocarbon', fontsize = 14)
		# y-axis title
		ax0.set_ylabel('H:C ratio of hydrocarbon', fontsize = 14)
		# set label font size
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

		# include colorbar to demonstrate time through experiment (s)
		cb = fig.colorbar(cm.ScalarMappable(norm=norm1, cmap=colors), ax=ax0)
		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('Normalised number concentration at this time', size=14, rotation=270, labelpad=20)


	return()

def plotter_mass_defect(self): # for mass defect plots

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = self.ro_obj.yrec
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = self.ro_obj.Nrec_wet
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = self.ro_obj.rad
	space_mode = self.ro_obj.spacing
	
	# get exact molar masses (g/mol)
	y_MW = np.array((y_MW))

	# get nominal mass (g/mol)
	nom_mass = np.array((ro_obj.nominal_mass))

	# calculate the mass defect according to doi.org/10.1002/jms.2953226 (g/mol)
	mass_def = y_MW-nom_mass

	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))

	if (self.oandm == 8): # all components in chemical scheme without scaling

		# plot mass defect against nominal mass
		ax0.plot(nom_mass, mass_def, '+')

	# all components in chemical scheme scaled to concentration at provided time
	if (self.oandm == 9):

		self.HC = np.array((ro_obj.HyC)) # hydrogen:carbon ratios of each component

		# get the time through experiment (s) in question
		try:
			treq = float((self.e403.toPlainText()))
		except:
			treq = timehr[-1]*3600.

		# index of the required experiment period
		tindx = np.min(np.abs(timehr*3600-treq)) == (np.abs(timehr*3600-treq))

		# index of hydrocarbons
		HCindx = self.HC > 0.

		# isolate gas-phase concentrations of hydrocarbons at this time (ppb)
		yrel = (yrec[tindx, 0:num_comp])[0, HCindx]

		# normalise to sum of concentrations
		yrel = yrel/sum(yrel)

		# colormap from library (https://matplotlib.org/stable/api/cm_api.html#matplotlib.cm.get_cmap)
		colors = cm.get_cmap('ocean')

		# set contour levels
		levels = (MaxNLocator(nbins = 256).tick_values(np.min(yrel), np.max(yrel)))
	
		# associate colours and contour levels
		norm1 = BoundaryNorm(levels, ncolors=256, clip=True)
	
		for i in range(len(yrel)):
			ax0.plot(nom_mass[HCindx][i], mass_def[HCindx][i], 'o', mec = 'k', mfc = colors(yrel[i]))

		# include colorbar to demonstrate time through experiment (s)
		cb = fig.colorbar(cm.ScalarMappable(norm=norm1, cmap=colors), ax=ax0)
		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('Normalised number concentration at this time', size=14, rotation=270, labelpad=20)

		ax0.set_title(str('Mass Defect at ' + str(treq) + ' s since experiment start for all hydrocarbons'), fontsize = 14)
			
	# x-axis title
	ax0.set_xlabel('Nominal Molar Mass (g/mol)', fontsize = 14)
	# y-axis title
	ax0.set_ylabel('Mass Defect (Exact Mass-Nominal Mass) (g/mol)', fontsize = 14)
	# set label font size
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	return()
	
def plotter_exp_prep(self): # for experiment design

	# get required variables from self
	wall_on = self.ro_obj.wf
	yrec = self.ro_obj.yrec
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	Nwet = self.ro_obj.Nrec_wet
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = self.ro_obj.rad
	space_mode = self.ro_obj.spacing
	
	# simulation title ---------------------------------------------
	sim_title = self.dir_path
	# -----------------------------------------------------------------
	
	# -----------------------------------------------------------------
	rationale = [] # remember to leave column for rationale
	# -----------------------------------------------------------------
	
	# -----------------------------------------------------------------
	# ratio of HO2, NO and RO2 reactivity with PINALO2 
	# averaged over whole simulation
	
	# prepare to hold values
	ct_PINALO2 = [0]*3
	
	# name of file that would contain change tendencies of PINALO2
	fname = str(self.dir_path+ '/PINALO2_rate_of_change')
	
	try: # try to open
		dydt = np.loadtxt(fname, delimiter = ',', skiprows = 0) # skiprows = 0 includes header
		# isolate chemical reaction numbers, note last two columns for gas-particle and gas-wall partitioning
		dydt_header = dydt[0, 0:-2]	
		dydt = dydt[1::, :] # exclude header now
		
	except: # if unable to open file
		PINALO2_ct = str('No change tendency record for PINALO2 was found, was it specified in the tracked_comp input of the model variables file?  Please see README for more information.')
	
	# prepare to store results of change tendency due to chemical reactions
	res = np.zeros((dydt.shape[0], dydt.shape[1]-2))
	res[:, :] = dydt[:, 0:-2] # get chemical reaction numbers and change tendencies
		
	import scipy.constants as si
	
	# get all reactions out of the used chemical scheme --------------------------------------------------------------------
	import sch_interr # for interpeting chemical scheme
	import re # for parsing chemical scheme
	import scipy.constants as si

	sch_name = self.ro_obj.sp
	inname = self.ro_obj.vp
	
	f_open_eqn = open(sch_name, mode= 'r' ) # open model variables file
	# read the file and store everything into a list
	total_list_eqn = f_open_eqn.readlines()
	f_open_eqn.close() # close file
		
	inputs = open(inname, mode= 'r' ) # open model variables file
	in_list = inputs.readlines() # read file and store everything into a list
	inputs.close() # close file
	
	# default chemical scheme markers
	self.chem_sch_mrk = ['{', 'RO2', '+', 'C(ind_', ')','' , '&', '' , '', ':', '}', ';']
	
	for i in range(len(in_list)): # loop through supplied model variables to interpret

		# ----------------------------------------------------
		# if commented out continue to next line
		if (in_list[i][0] == '#'):
			continue
		key, value = in_list[i].split('=') # split values from keys
		# model variable name - a string with bounding white space removed
		key = key.strip()
		# ----------------------------------------------------

		if key == 'chem_scheme_markers' and (value.strip()): # formatting for chemical scheme
			self.chem_sch_mrk = [str(i).strip() for i in (value.split(','))]

	# interrogate scheme to list equations
	[eqn_list, aqeqn_list, eqn_num, rrc, rrc_name, 
		RO2_names] = sch_interr.sch_interr(total_list_eqn, self)	
	
	# list interested reactants
	reac_interest = [' HO2 ', '*RO2', ' NO ']
	
	eqn_cnt = 0 # count on equations through chemical scheme
	PINALO2_eqn_cnt = 0 # count on PINALO2 equations
	
	# loop through equations to find where PINALO2 reacts with HO2, RO2 or NO
	for eqn_numi in dydt_header:
		
		# obtain original equation
		eqni = eqn_list[int(eqn_numi)]
		
		for reaci in range(3): # loop through the reactants we're interested in
		
			# focus on LHS of reaction first
			if reac_interest[reaci] in eqni.split('=')[0]:
				# sum change tendency over time and add to total
				ct_PINALO2[reaci] += np.sum(dydt[:, PINALO2_eqn_cnt])
				break # move onto next equation
			
			# the reaction rate coefficient for this equation
			# where reaction rate coefficient starts
			indx_rrc_start = eqni.index(self.chem_sch_mrk[9])
			# where equation starts
			indx_rrc_end = eqni.index(self.chem_sch_mrk[10])
			
			# in case equation comes before the reaction rate coefficient,
			# then change end point of reaction rate coefficient part to end 
			# of line
			if (indx_rrc_end < indx_rrc_start):
				indx_rrc_end = eqni.index(self.chem_sch_mrk[11])
				
			# likewise, if reactant is in the reaction rate coefficient
			if reac_interest[reaci] in eqni[indx_rrc_start:indx_rrc_end]:
				
				# sum change tendency over time and add to total
				ct_PINALO2[reaci] += np.sum(dydt[:, PINALO2_eqn_cnt])
				break
		
		PINALO2_eqn_cnt+= 1 # count on equations

	# -----------------------------------------------------------------

	# open an excel file ready for saving results
	from openpyxl import Workbook
	wb = Workbook() # create new workbook
	ws = wb.active # get active worksheet
	ws.title = "Exp. Prep." # set worksheet title
	# add headers
	ws["A1"] = 'Exp. #'
	ws["B1"] = 'Rationale'
	ws["C1"] = 'HO2:NO:RO2 reactivity with PINALO2 over whole simulation'
	# add data
	ws["A2"] = 1
	ws["B2"] = 'HO2 and RO2 at play'
	ws["C2"] = str(str(ct_PINALO2[0]) + ':' + str(ct_PINALO2[2]) + ':' + str(ct_PINALO2[1]))
	wb.save(filename = str(self.dir_path+ '/exp_prep.xlsx')) # save workbook
	
	
	return() 