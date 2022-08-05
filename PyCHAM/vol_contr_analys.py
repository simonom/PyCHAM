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
'''script for plotting volatility basis set mass fraction of particle phase with time and tabulating component volatilities and particle-phase concentrations'''
# aids interpretation of gas-particle partitioning results, 
# assumes calling from the PyCHAM home directory

# import required modules
import os
import sys
# ensure modules can be seen 
# (assumes calling from the home folder)
sys.path.append(str(os.getcwd() + '/PyCHAM'))
import numpy as np
import pybel
import xml_interr
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap # for customised colormap
import matplotlib.ticker as ticker # set colormap tick labels to standard notation
import scipy.constants as si
import retr_out

def plotter_wiw(caller, dir_path, self, now): # define function

	# inputs: -------------------------------
	# caller - the module calling (0 for gui)
	# dir_path - path to results
	# self - reference to PyCHAM
	# now - whether to include (0) or exclude water (1)
	# self.phase4vol - the phase to consider for plotting results
	# -----------------------------------------

	# prepare the volatility basis set interpretation
	# of particle-phase concentrations

	# required outputs from full-moving
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, rel_SMILES, 
		y_mw, N, comp_names, y_MV, _, wall_on, space_mode, _, _, _, 
		PsatPa, OC, H2Oi, seedi, _, _, _, _, _) = retr_out.retr_out(dir_path)
	

	# tell user that this code won't work if no particle size bins present
	if (self.phase4vol == 'Particle Phase' or self.phase4vol == 'Particle and Gas Phase Combined'):
		
		if ((num_sb-wall_on) == 0 or sum(sum(Ndry))== 0.):
			self.l203a.setText('Error - volatility basis set for particle phase may only be estimated for simulations where gas-particle partitioning has occurred')
			if (self.bd_pl == 1):
				self.l203a.setStyleSheet(0., '2px dashed red', 0., 0.)
				self.bd_pl = 2
			else:
				self.l203a.setStyleSheet(0., '2px solid red', 0., 0.)
				self.bd_pl = 1
			return()
		
	# prepare plot -----------------------------------------------------------------------
	if (caller == 0): # if calling function is gui
		plt.ion() # show figure
	
	# prepare plot
	fig, (ax1) = plt.subplots(1, 1, figsize=(10,7))
	fig.subplots_adjust(hspace = 0.7)

	# ----------------------------------------------------------------------------------------

	# number of particle size bins without wall
	num_asb = (num_sb-wall_on)
	
	# convert from list to array
	y_mw = (np.array((y_mw))).reshape(1, -1)
	PsatPa = (np.array((PsatPa))).reshape(1, -1)
	
	# repeat molecular weights over size bins and times
	y_mw_rep = np.tile(y_mw, (1, num_asb))
	y_mw_rep = np.tile(y_mw_rep, (len(t_array), 1))
	
	if (self.phase4vol == 'Particle Phase'):
	
		# particulate concentrations of individual components (*1.e-12 to convert from g/cm3 (air) to ug/m3 (air))
		# including any water and core
		pc = (y[:, num_comp:num_comp*(num_asb+1)]/si.N_A)*y_mw_rep*1.e12
	
		if (now == 1): # if water to be excluded, zero its contribution
			pc[:, H2Oi::num_comp] = 0.
	
	if (self.phase4vol == 'Gas Phase'):
		# gaseous concentrations of individual components (*1.e-12 to convert from g/cm3 (air) to ug/m3 (air))
		# including any water and core
		pc = (y[:, 0:num_comp]*Cfac/si.N_A)*y_mw*1.e12
		num_asb = 1
		
		if (now == 1): # if water to be excluded, zero its contribution
			pc[:, H2Oi::num_comp] = 0.
		
	if (self.phase4vol == 'Gas Phase Only C>1, O>0'):
		# gaseous concentrations of individual hydrocarbon components (*1.e-12 to convert from g/cm3 (air) to ug/m3 (air))
		# including any water and core
		pc = (y[:, 0:num_comp]*np.array((Cfac)).reshape(-1, 1)/si.N_A)*np.array((y_mw)).reshape(1, -1)*1.e12
		
		# prepare for results
		hc_indx = []
		
		for smile_now in rel_SMILES:
			hc_indx.append(smile_now.count('C')>1*smile_now.count('O')>0)
		
		pc = pc[:, hc_indx]
		PsatPa = PsatPa[0, hc_indx]
		
		y_mw = y_mw[0, hc_indx]

		num_asb = 1

	
	
	# total particulate concentrations at each time (ug/m3)
	tpc = pc.sum(axis=1)
	
	# standard temperature for pure component saturation vapour pressures (K)
	TEMP = 298.15

	# convert standard (at 298.15 K) vapour pressures in Pa to 
	# saturation concentrations in ug/m3
	# using eq. 1 of O'Meara et al. (2014)
	Psat_Cst = (1.e6*y_mw[0, :])*(PsatPa/101325.)/(8.2057e-5*TEMP)
	
	# tile over size bins
	Psat_Cst  = np.tile(Psat_Cst, num_asb)
	# remove excess dimension
	Psat_Cst = Psat_Cst.squeeze()

	# the saturation concentrations to consider (log10(C* (ug/m3)))
	# note these will be values at the centre of the volatility size bins
	# setting the final input argument to 1 means decadal bins of vapour pressure
	sc = np.arange(-2.5, 7.5, 1.)

	# empty array for normalised mass contributions
	nmc = np.zeros((len(sc), len(t_array)))

	for it in range(len(t_array)): # loop through times
	
		# loop through saturation concentrations and find normalised mass contributions
		# to particulate loading
		for i in range(len(sc)):

			if (i == 0):
				indx = Psat_Cst < 10**(sc[i]+0.5)
			if (i > 0 and i < len(sc)-1):
				indx = (Psat_Cst >= 10**(sc[i-1]+0.5))*(Psat_Cst < 10**(sc[i]+0.5))
			if (i == len(sc)-1):
				indx = (Psat_Cst >= 10**(sc[i-1]+0.5))

			if (tpc[it] > 0.):
				nmc[i, it] = (pc[it, indx].sum())/tpc[it]

	
	# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
	colors = [(0.6, 0., 0.7), (0, 0, 1), (0, 1., 1.), (0, 1., 0.), (1., 1., 0.), (1., 0., 0.)]  # R -> G -> B
	n_bin = 100  # discretizes the colormap interpolation into bins
	cmap_name = 'my_list'
	# create the colormap
	cm = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
	
	# set contour levels
	levels = (MaxNLocator(nbins = 100).tick_values(np.min(nmc), np.max(nmc)))
	
	# associate colours and contour levels
	norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)

	ptindx = (tpc > 0) # indices of times where secondary material present

	p0 = ax1.pcolormesh(t_array[ptindx], sc, nmc[:, ptindx], cmap=cm, norm=norm1, shading='auto')

	cax = plt.axes([0.875, 0.40, 0.02, 0.18]) # specify colour bar position
	cb = plt.colorbar(p0, cax = cax, ticks=[0.00, 0.25, 0.50, 0.75, 1.00], orientation = 'vertical')
	cb.ax.tick_params(labelsize = 12)
	cb.set_label('mass fraction', size = 12, rotation = 270, labelpad = 10.)

	ax1.set_xlabel(r'Time through experiment (hours)', fontsize=14)
	ax1.set_ylabel(r'$\rm{log_{10}(}$$C*_{\mathrm{298.15 K}}$$\rm{\, (\mu g\, m^{-3}))}$', fontsize=14, labelpad = 10.)
	ax1.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax1.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	# array containing the location of tick labels
	ytloc = sc
	ax1.set_yticks(ytloc)
	# prepare list of strings for the tick labels
	ytl = []
	for i in ytloc:
		if (i == np.min(ytloc)): # if the minimum include less than sign
			ytl.append(str('$\less$' +str(i+0.5)))
			continue
		if (i == np.max(ytloc)): # if the maximum include the greater than or equal to sign
			ytl.append(str('$\geq$' +str(i+0.5)))
			continue
		ytl.append(str(i+0.5)) # otherwise just state number
		
	ax1.set_yticklabels(ytl)

	if (caller != 0):
		plt.show() # show figure

	return() # end function

# define function for plotting the two-dimensional volatility basis set,
# e.g. https://doi.org/10.5194/acp-20-1183-2020
def plotter_2DVBS(caller, dir_path, self, t_thro):

	# inputs: -------------------------------
	# caller - the module calling (0 for gui)
	# dir_path - path to results
	# self - reference to GUI
	# t_thro - time (s) through experiment at which to pot the 2D-VBS
	# -----------------------------------------

	# ----------------------------------------------------------------------------------------
	if (caller == 0): # if calling function is gui
		plt.ion() # show figure
	
	# prepare plot
	fig, (ax1) = plt.subplots(1, 1, figsize=(10,7))
	fig.subplots_adjust(hspace = 0.7)
	
	# prepare plot data --------------------------------------
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, rel_SMILES, 
		y_mw, N, comp_names, y_MV, _, wall_on, space_mode, _, _, _, PsatPa, OC, 
		H2Oi, seedi, _, _, _, _, _) = retr_out.retr_out(dir_path)
	
	# subtract recorded times from requested time and absolute
	t_diff = np.abs(t_thro-(t_array*3600.))
	
	# find closest recorded time to requested time to plot
	t_indx = (np.where(t_diff == np.min(t_diff)))[0][0]
	
	# temperature for vapour pressures
	TEMP = 298.15
	
	# convert lists to numpy array
	y_mw = np.array((y_mw))
	PsatPa = np.array((PsatPa))
	
	# convert standard (at 298.15 K) vapour pressures in Pa to 
	# saturation concentrations in ug/m3
	# using eq. 1 of O'Meara et al. 2014
	Psat_Cst = (1.e6*y_mw)*(PsatPa/101325.)/(8.2057e-5*TEMP)
	
	# get particle concentrations at this time (# molecules/cm3)
	pc = y[t_indx, num_comp:num_comp*(num_sb-wall_on)] 

	# zero water and seed components
	pc[H2Oi[0]::num_comp] = 0.
	for seed_indxs in seedi:
		pc[seed_indxs::num_comp] = 0.
	
	# tile molecular weights over size bins
	y_mw = np.tile(y_mw, (num_sb-wall_on-1))
	
	# convert concentrations from # molecules/cm3 to ug/m3
	pc = ((pc/si.N_A)*y_mw)*1.e12
	
	# sum particle concentrations (ug/m3) at this time, without water and seed
	tot_pc = pc.sum()
	
	OC_range = np.arange(0., 2., 0.2) # the O:C range
	VP_range = np.arange(-2.5, 7.5, 1.) # the log10 of the vapour pressure (ug/m3) at 298.15 K range
	
	# empty mass fractions matrix
	mf = np.zeros((len(OC_range), len(VP_range)))
	
	# convert list to array
	PsatPa = np.array((PsatPa))
	OC = np.array((OC))
	
	# loop over size bins to get total particle-phase 
	# concentration of each component (ug/m3)
	for sbi in range(1, (num_sb-wall_on)-1):
		pc[0:num_comp] += pc[num_comp*sbi:num_comp*(sbi+1)]
	
	# forget all excess size bins (ug/m3)
	pc = pc[0:num_comp]
	
	# loop through O:C ratios and vapour pressures to estimate mass fractions
	for OCi in range(len(OC_range)-1):
	
		VPo = 0. # reset lower vapour pressure limit (ug/m3)
		
		for VPi in range(len(VP_range)):
			
			# upper vapour pressure limit now (ug/m3)
			VPn = 10**(VP_range[VPi]+0.5)
			if (VPi == len(VP_range)-1): # final upper vapour pressure limit (ug/m3)
				VPn = np.inf
			if (OCi == len(OC_range)-1): # final upper O:C
				OC_up = np.inf
			else:
				OC_up = OC_range[OCi+1]
				
			# index of all components with this vapour pressure and O:C ratio
			compi =  ((Psat_Cst >= VPo)*(Psat_Cst < VPn))*((OC >= OC_range[OCi])*(OC < OC_up))

			# mass fractions of components with this combination of properties
			if (tot_pc > 0):
				mf[OCi, VPi] = sum(pc[compi])/tot_pc

			VPo = VPn # reset lower vapour pressure limit (ug/m3)
	
	# do the plotting ------------------------
	
	# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
	colors = [(0.6, 0., 0.7), (0, 0, 1), (0, 1., 1.), (0, 1., 0.), (1., 1., 0.), (1., 0., 0.)]  # R -> G -> B
	n_bin = 100  # discretizes the colormap interpolation into bins
	cmap_name = 'my_list'
	# create the colormap
	cm = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
	
	# set contour levels
	levels = (MaxNLocator(nbins = 100).tick_values(0., 1.))
	
	# associate colours and contour levels
	norm1 = BoundaryNorm(levels, ncolors = cm.N, clip=True)

	p0 = ax1.pcolormesh(VP_range, OC_range, mf, cmap=cm, norm=norm1, shading='auto')

	cax = plt.axes([0.875, 0.40, 0.02, 0.18]) # specify colour bar position
	cb = plt.colorbar(p0, cax = cax, ticks=[0.00, 0.25, 0.50, 0.75, 1.00], orientation = 'vertical')
	cb.ax.tick_params(labelsize = 12)
	cb.set_label('mass fraction', size = 12, rotation = 270, labelpad = 10.)

	ax1.set_xlabel(r'$\rm{log_{10}(}$$C*_{\mathrm{298.15 K}}$$\rm{\, (\mu g\, m^{-3}))}$', fontsize=14)
	ax1.set_ylabel(r'O:C ratio', fontsize=14, labelpad = 10.)
	ax1.set_title(str('Mass fraction of non-water and non-seed components at ' + str(t_thro) + str(' s through experiment')), fontsize=14)
	
	ax1.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax1.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	
	# array containing the location of tick labels
	xtloc = VP_range
	ax1.set_xticks(xtloc)
	ytloc = OC_range
	ax1.set_yticks(ytloc)
	# prepare list of strings for the tick labels
	xtl = []
	for i in xtloc:
		if (i == np.min(xtloc)): # if the minimum include less than sign
			xtl.append(str('$\less$' +str(i+0.5)))
			continue
		if (i == np.max(xtloc)): # if the maximum include the greater than or equal to sign
			xtl.append(str('$\geq$' +str(i+0.5)))
			continue
		xtl.append(str(i+0.5)) # otherwise just state number
		
	ax1.set_xticklabels(xtl)

	if (caller != 0):
		plt.show() # show figure

	return() # end function
	
# define function for plotting gas-phase concentrations at the 
# user-defined time through experiment ordered by volatility
def plotter_gpc(caller, dir_path, self, t_thro):

	# inputs: -------------------------------
	# caller - the module calling (0 for gui)
	# dir_path - path to results
	# self - reference to GUI
	# t_thro - time (s) through experiment at which to pot the 2D-VBS
	# -----------------------------------------

	# ----------------------------------------------------------------------------------------
	if (caller == 0): # if calling function is gui
		plt.ion() # show figure
	
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(10,7))
	fig.subplots_adjust(hspace = 0.7)
	
	# prepare plot data --------------------------------------
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, rel_SMILES, 
		y_mw, N, comp_names, y_MV, _, wall_on, space_mode, _, _, _, PsatPa, OC, 
		H2Oi, seedi, _, _, _, _, _) = retr_out.retr_out(dir_path)
	
	# subtract recorded times from requested time and absolute
	t_diff = np.abs(t_thro-(t_array*3600.))
	
	# find closest recorded time to requested time to plot
	t_indx = (np.where(t_diff == np.min(t_diff)))[0][0]
	
	# isolate gas-phase concentrations (ppb) at this time
	y_gp_ppb = np.squeeze(y[t_indx, 0:num_comp])

	# order components by their pure component saturation vapour pressures at 298.15 K
	# get indices of vapour pressure in descending order
	des_ind = np.flip(np.argsort(PsatPa, axis = 0))
	
	array_names = np.squeeze(np.array(comp_names))[des_ind]
	VPs = np.array((PsatPa))[des_ind]
	
	# component names and their vapour pressures listed in descending order of volatility 
	# together and with 1 decimal point precision
	xticks_str = [] # start with empty list
	
	for ci in range(len(array_names)):
		xticks_str.append(str(array_names[ci] + ' (' + "{:.2e}".format(VPs[ci]) + ')'))
	
	
	# plot gas-phase concentrations against component names in descending order of volatility
	ax0.semilogy(np.arange(len(y_gp_ppb)), y_gp_ppb[des_ind], '+')

	ax0.set_ylabel(str('Gas-phase Concentration at ' + str(t_thro) + ' s Through Experiment (ppb)'), fontsize = 14)
	ax0.set_xlabel(r'Component Name (Pure Component Saturation Vapour Pressure at 298.15 K (Pa))', fontsize = 14)
	# set location of x ticks
	ax0.set_xticks(np.arange(len(y_gp_ppb)))
	ax0.set_xticklabels(xticks_str, rotation = 45)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')

	return()
	
def plotter_pie_top_n(self): # define function to plot the top n mass contributors to a phase at a given time

	# inputs: -------------------------------
	# self - reference to GUI
	# -----------------------------------------

	# ----------------------------------------------------------------------------------------
	plt.ion() # show figure
	
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(10,7))
	fig.subplots_adjust(hspace = 0.7)
	
	# prepare plot data --------------------------------------
	# required outputs from full-moving
	(num_sb, num_comp, Cfac, y, Ndry, rbou_rec, xfm, t_array, rel_SMILES, 
		y_mw, N, comp_names, y_MV, _, wall_on, space_mode, _, _, _, PsatPa, OC, 
		H2Oi, seedi, _, _, _, _, _) = retr_out.retr_out(self.dir_path)
	
	# subtract recorded times from requested time and absolute
	t_diff = np.abs(self.t_thro-(t_array*3600.))
	
	# find closest recorded time to requested time to plot
	t_indx = (np.where(t_diff == np.min(t_diff)))[0][0]
	
	comp_names = np.array((comp_names)) # ensure array
	PsatPa = np.array((PsatPa)) # ensure array

	y_mw = np.array((y_mw)) # ensure array

	# particle-phase
	# isolate concentrations for phase to consider
	if (self.phase4vol == 'Particle Phase'):
		import scipy.constants as si
		yn = y[t_indx, num_comp:(num_sb-wall_on+1)*num_comp] # particle-phase concentrations (# molecules/cm3)
		# convert # molecules/cm3 to ug/m3
		yn = ((yn)/(si.N_A))*y_mw*1.e12

		# sum over size bins
		yn = np.squeeze(yn.reshape((num_sb-wall_on), num_comp).sum(axis=0))

	# gas-phase, for non-methane oxidised components
	if (self.phase4vol == 'Gas Phase Only C>1, O>0'):
		import scipy.constants as si
		yn = y[t_indx, 0:num_comp] # gas-phase concentrations (ppb)
		# convert ppb to ug/m3
		yn = ((yn*Cfac[t_indx])/(si.N_A))*y_mw*1.e12
	
		# now extract only components with more than one carbon and containing oxygen
		# prepare for results
		hc_indx = []
		
		for smile_now in rel_SMILES:
			hc_indx.append(smile_now.count('C')>1*smile_now.count('O')>0)
			
		yn = yn[hc_indx]

		comp_names = comp_names[hc_indx]
		
		PsatPa = PsatPa[hc_indx]

		y_mw = y_mw[hc_indx]

	# get mass fractions
	yn = yn/sum(yn)
	
	# ascending order by mass fractions
	asc_ind = np.argsort(yn, axis = 0)
	
	# get top number of fractions
	yn = yn[asc_ind][-self.num_pie_comp::]

	# include fraction composed of other components
	yn = np.concatenate((yn, np.array((1.-sum(yn)).reshape(1))))
	
	# get top component names
	top_name = comp_names[asc_ind][-self.num_pie_comp::]
	
	# get top component vapour pressures (Pa)
	top_vp = PsatPa[asc_ind][-self.num_pie_comp::]

	# get top component molar masses (g/mol)
	top_mw = y_mw[asc_ind][-self.num_pie_comp::]

	# prepare for pie section labels
	label = []
	
	# convert standard (at 298.15 K) vapour pressures in Pa to 
	# saturation concentrations in ug/m3 (at 298.15 K)
	# using eq. 1 of O'Meara et al. (2014)
	Psat_Cst = (1.e6*top_mw)*(top_vp/101325.)/(8.2057e-5*298.15)

	# combine names and vapour pressures
	for i in range(len(top_name)):
		label.append((str(top_name[i] + ' (' + "{:.2e}".format(Psat_Cst[i]) + r' $\rm{\mu g\, m^{-3})}$')))
	
	# final label is for all other components
	label.append('remainder')

	# plot pie chart with labels and title
	ax0.pie(yn, labels = label, normalize=True)
	
	# include title
	ax0.set_title(str('Contribution to mass of ' + self.phase4vol + ' at ' + str(self.t_thro) + ' s through experiment for top ' + str(self.num_pie_comp) + ' components (C*)' ), fontsize=14)
	
	return()