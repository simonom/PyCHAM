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
'''standard graphical representation of simulation results'''
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

def plotter(caller, dir_path, uc, self):
	
	# inputs: ------------------------------------------------------------------
	# caller - marker for whether PyCHAM (0) or tests (2) are the calling module
	# dir_path - path to folder containing results files to plot
	# uc - number representing the units to be used for gas-phase concentrations
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
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing

	# number of actual particle size bins
	num_asb = (num_sb-wall_on)

	if (caller == 0):
		plt.ion() # show results to screen and turn on interactive mode

	# prepare sub-plots depending on whether particles present
	if (num_asb == 0): # no particle size bins
		if not (indx_plot): # check whether there are any gaseous components to plot
			mess = str('Please note, no initial gas-phase concentrations were received and no particle size bins were present, therefore there is nothing for the standard plot to show')
			self.l203a.setText(mess)
			return()

		# if there are gaseous components to plot, then prepare figure
		fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
		mess = str('Please note, no particle size bins were present, therefore the particle-phase standard plot will not be shown')
		self.l203a.setText(mess)

	if (num_asb > 0): 
		
		if not (indx_plot): # no gaseous components
			mess = str('Please note, no initial gas-phase concentrations were registered, therefore the gas-phase standard plot is not shown')
			self.l203a.setText(mess)

			
		if not (indx_plot): # no gaseous components, only particle-phase
			mess = str('Please note, no initial gas-phase concentrations were registered, therefore the gas-phase standard plot will not be shown')
			self.l203a.setText(mess)
			# if there are no gaseous components then prepare figure
			fig, (ax1) = plt.subplots(1, 1, figsize=(14, 7))
			
				
		if (indx_plot):
			# if there are both gaseous components and particle 
			# size bins then prepare figure
			fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 7))

		# parasite axis setup on particle-phase plot-----------------------------
		par1 = ax1.twinx() # first parasite axis
		par2 = ax1.twinx() # second parasite axis
		# Offset the right spine of par2.  The ticks and label have already been
		# placed on the right by twinx above.
		par2.spines["right"].set_position(("axes", 1.2))
		# Having been created by twinx, par2 has its frame off, so the line of its
		# detached spine is invisible.  First, activate the frame but 
		# make the patch and spines invisible.
		make_patch_spines_invisible(par2)
		# second, show the right spine
		par2.spines["right"].set_visible(True)	
		# -----------------------------------------------------------------------

	if (indx_plot):
		
		
		ymax = 0. # start tracking maximum value for plot label
		
		# action units for gas-phase concentrations
		if (uc == 0): # ppb
			gp_conc = yrec[:, 0:num_comp] # ppb is original units
			gpunit = '(ppb)'
		if (uc == 1 or uc == 2): # ug/m3 or # molecules/cm3

			y_MW = (np.array(y_MW)).reshape(1, -1) # convert to numpy array from list
			Cfaca = (np.array(Cfac)).reshape(-1, 1) # convert to numpy array from list
			
			gp_conc = yrec[:, 0:num_comp] 

			# # molecules/cm3
			gp_conc = gp_conc.reshape(yrec.shape[0], num_comp)*Cfaca
			gpunit = str('\n(' + u'\u0023' + ' molecules/cm' + u'\u00B3' + ')')

			if (uc == 1): # ug/m3
				gp_conc = ((gp_conc/si.N_A)*y_MW)*1.e12
				gpunit = str('(' + u'\u03BC' + 'g/m' +u'\u00B3' + ')')

		# gas-phase concentration sub-plot ---------------------------------------------	
		for i in range(len(indx_plot)):
		
			ax0.semilogy(timehr, gp_conc[:, indx_plot[i]], '+',linewidth=4.0, 
						label=str(str(comp0[i]).strip()))
			ymax = max(ymax, max(yrec[:, indx_plot[i]]))
		if (uc == 1 or uc == 2): # ug/m3 or # molecules/cm3
			ax0.set_ylabel(r'Gas-phase concentration ' + gpunit, fontsize = 14)
		if (uc == 0): # ppb
			ax0.set_ylabel(r'Gas-phase mixing ratio ' + gpunit, fontsize = 14)
		ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
		ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax0.legend(fontsize=14)

		# find maximum and minimum of plotted concentrations for sub-plot label		
		maxy = max(yrec[:, indx_plot].flatten())
		miny = min(yrec[:, indx_plot].flatten())	
		
		if (num_asb > 0): # if more than one plot
			# get the location of ticks
			locs = ax0.get_yticks()
			maxloc = max(locs)
			ax0.text(x = timehr[0]-(timehr[-1]-timehr[0])/9.5, 
				y = ymax*1.05, s='a)', size=14)

		# end of gas-phase concentration sub-plot -----------------------------------
	
	# particle properties sub-plot --------------------------------------------------
	if (num_asb > 0): # if size bins present		

		if (timehr.ndim == 0): # occurs if only one time step saved
			Nwet = np.array(Nwet.reshape(1, num_asb))
			rbou_rec = np.array(rbou_rec.reshape(1, num_sb))
		if (num_asb == 1): # just one particle size bin (wall included in num_sb)
			Nwet = np.array(Nwet.reshape(len(timehr), num_asb))

		# plotting number size distribution --------------------------------------
	
		# don't use the first boundary as it could be zero, which will error when log10 taken
		log10D = np.log10(rbou_rec[:, 1::]*2.)
		
		if (num_asb > 1) :
			# note, can't append zero to start of log10D to cover first size bin as the log10 of the
			# non-zero boundaries give negative results due to the value being below 1, so instead
			# assume same log10 distance as the next pair
			log10D = np.append((log10D[:, 0]-(log10D[:, 1]-log10D[:, 0])).reshape(-1, 1), log10D, axis=1)
			# radius distance covered by each size bin (log10(um))
			dlog10D = (log10D[:, 1::]-log10D[:, 0:-1]).reshape(log10D.shape[0], log10D.shape[1]-1)
		
		if (num_asb == 1): # single particle size bin
			# assume lower radius bound is ten times smaller than upper
			dlog10D = (log10D[:, 0]-np.log10((rbou_rec[:, 1]/10.)*2.)).reshape(log10D.shape[0], 1)
				
		# number size distribution contours (# particles/cm3 (air))
		dNdlog10D = np.zeros((Nwet.shape[0], Nwet.shape[1]))
		dNdlog10D[:, :] = Nwet[:, :]/dlog10D[:, :]
		# transpose ready for contour plot
		dNdlog10D = np.transpose(dNdlog10D)
		
		# mask any nan values so they are not plotted
		z = np.ma.masked_where(np.isnan(dNdlog10D), dNdlog10D)
	
		# customised colormap (https://www.rapidtables.com/web/color/RGB_Color.html)
		colors = [(0.6, 0., 0.7), (0, 0, 1), (0, 1., 1.), (0, 1., 0.), (1., 1., 0.), (1., 0., 0.)]  # R -> G -> B
		n_bin = 100  # discretizes the colormap interpolation into bins
		cmap_name = 'my_list'
		# create the colormap
		cm = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
	
		
		# -----------------------------
		# smallest value to consider for levels
		z_min = np.max(z[~np.isnan(z)])*1.e-3
	
		if np.max(z[~np.isnan(z)]) == 0: # if no particle present above 0 /cm3
			levels = np.zeros((1))
			
		else:	
			# make a first contour plot (which will be covered by plot p1 below) to get 
			# the contour labels of actual concentration (not log10(concentration))
			# set contour levels
			if (z_min > 0.5): # if rounding gives a number greater than 0 
				levels = np.arange(np.log10(round(z_min)), np.log10(np.max(z[~np.isnan(z)])), (np.log10(np.max(z[~np.isnan(z)]))-np.log10(round(z_min)))/1.e2)
			else: # i rounding would give a number below zero, then don't round		
				levels = np.arange(np.log10((z_min)), np.log10(np.max(z[~np.isnan(z)])), (np.log10(np.max(z[~np.isnan(z)]))-np.log10((z_min)))/1.e2)
			# associate colours and contour levels
			norm1 = BoundaryNorm(10.**levels, ncolors=cm.N, clip=True)
		
			# contour plot with times (hours) along x axis and 
			# particle diameters (nm) along y axis
			for ti in range(len(timehr)-1): # loop through times
				p0 = ax1.pcolormesh(timehr[ti:ti+2], (rbou_rec[ti, :]*2*1e3), (z[:, ti]).reshape(-1, 1), cmap=cm, norm=norm1)
		
		# ----------------------------
		
		if np.max(z[~np.isnan(z)]) == 0: # if no particle present above 0 /cm3
			
			levels = np.arange(-0.1, 0.1, (0.1--0.1)/1.e2)
			# associate colours and contour levels
			norm1 = BoundaryNorm(levels, ncolors=cm.N, clip=True)
			
			# contour plot with times (hours) along x axis and 
			# particle diameters (nm) along y axis
			for ti in range(len(timehr)-1): # loop through times
				p1 = ax1.pcolormesh(timehr[ti:ti+2], (rbou_rec[ti, :]*2*1e3), (z[:, ti]).reshape(-1, 1), cmap=cm, norm=norm1)

			cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25, ax=ax1)
			
		else:
			if (z_min > 0.5): # if rounding would give zero
				# set contour levels
				levels = np.arange(np.log10(round(z_min)), np.log10(np.max(z[~np.isnan(z)])), (np.log10(np.max(z[~np.isnan(z)]))-np.log10(round(z_min)))/1.e2)
			else: # don't round if minimum close to zero
				levels = np.arange(np.log10((z_min)), np.log10(np.max(z[~np.isnan(z)])), (np.log10(np.max(z[~np.isnan(z)]))-np.log10((z_min)))/1.e2)
			
			# associate colours and contour levels
			norm1 = BoundaryNorm(levels, ncolors=cm.N, clip=True)
		
			# get indices of zeros in z
			zero_indx = z == 0.
			# get minimum value above 0 in z
			z_gt_zero = np.min(z[z!=0])
			# temporarily assign a > 0 number to zeros
			z[zero_indx] = z_gt_zero*1.e-1
			z_log10 = np.log10(z)
			z[zero_indx] = 0. # return to zero
		
			# contour plot with times (hours) along x axis and 
			# particle diameters (nm) along y axis
			for ti in range(len(timehr)-1): # loop through times
				p1 = ax1.pcolormesh(timehr[ti:ti+2], (rbou_rec[ti, :]*2*1e3), (z_log10[:, ti]).reshape(-1, 1), cmap=cm, norm=norm1)
	
			cb = plt.colorbar(p0, format=ticker.FuncFormatter(fmt), pad=0.25, ax=ax1)
	
		# if logarithmic or manual spacing of size bins specified, plot vertical axis 
		# logarithmically
		if space_mode == 'log' or space_mode == 'man':
			ax1.set_yscale("log")
			
		# set tick format for vertical axis
		ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e'))
		ax1.set_ylabel('Diameter (nm)', size = 14)
		ax1.xaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')
		ax1.yaxis.set_tick_params(labelsize = 14, direction = 'in', which = 'both')

		# label according to whether gas-phase plot also displayed		
		if (indx_plot):
			ax1.text(x=timehr[0]-(timehr[-1]-timehr[0])/11., y = np.amax(rbou_rec*2*1e3)*1.05, s='b)', size=14)
		ax1.set_xlabel(r'Time through simulation (hours)', fontsize=14)
		
		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('dN (#$\,$$\mathrm{cm^{-3}}$)/d$\,$log$_{10}$(D$\mathrm{_p}$ ($\mathrm{\mu m}$))', size=14, rotation=270, labelpad=20)

		# ------------------------------------------------------------------
		# total particle number concentration # particles/cm3
	
		# include total number concentration (# particles/cm3 (air)) on contour plot
		# first identify size bins with radius exceeding 3nm
		# empty array for holding total number of particles
		Nvs_time = np.zeros((Nwet.shape[0]))
	
		for i in range(num_asb): # size bin loop
			Nvs_time[:] += Nwet[:, i] # sum number
		
		p3, = par1.plot(timehr, Nvs_time, '+k', label = 'N')
	
		par1.set_ylabel('N (#$\,$ $\mathrm{cm^{-3})}$', size=14, rotation=270, labelpad=20) # vertical axis label
		par1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e')) # set tick format for vertical axis
		par1.yaxis.set_tick_params(labelsize=14)

		# mass concentration of particles --------------------------------------------
		# array for mass concentration with time
		MCvst = np.zeros((1, len(timehr)))
		
		# first obtain just the particle-phase concentrations (# molecules/cm3)
		yrp = np.zeros((yrec.shape[0], num_comp*(num_asb)))
		yrp[:, :] = yrec[:, num_comp:num_comp*(num_asb+1)]
		# loop through size bins to convert to ug/m3
		for sbi in range(num_asb):
			yrp[:, sbi*num_comp:(sbi+1)*num_comp] = ((yrp[:, sbi*num_comp:(sbi+1)*num_comp]/si.N_A)*y_MW)*1.e12
		
		MCvst[0, :] = yrp.sum(axis=1)

		# log10 of maximum in mass concentration
		if (max(MCvst[0, :]) > 0):
			MCmax = int(np.log10(max(MCvst[0, :])))
		else:
			MCmax = 0.
	
		p5, = par2.plot(timehr, MCvst[0, :], 'xk', label = 'Total Particle Mass Concentration')
		par2.set_ylabel(str('Mass Concentration ($\mathrm{\mu g\, m^{-3}})$'), rotation=270, size=16, labelpad=25)
		# set colour of label, tick font and corresponding vertical axis to match scatter plot presentation
		par2.yaxis.label.set_color('black')
		par2.tick_params(axis='y', colors='black')
		par2.spines['right'].set_color('black')
		par2.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.1e')) # set tick format for vertical axis
		par2.yaxis.set_tick_params(labelsize=16)
		plt.legend(fontsize=14, handles=[p3, p5] , loc=4, fancybox=True, framealpha=0.5)	

	# end of particle properties sub-plot -----------------------------------

	if (caller == 2): # display when in test mode
		plt.show()	

	return()
	

def aqi_calc(self): # define function
	
	# inputs ----------------------------------------------------------
	# self - reference to PyCHAM class
	# -------------------------------------------------------------------
	
	# import dependencies
	import math
	import numpy as np
	import scipy.constants as si
	
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
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec= np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	Cfac = self.ro_obj.cfac
	cen_size = self.ro_obj.cen_size

	# time in minutes
	timemin = timehr*60.

	# limit of index for 24-hour running means
	indx_lim24 = sum(timehr <= (timehr[-1]-24.))

	# limit of index for 15-minute running means
	indx_lim15 = sum(timemin <= (timemin[-1]-15.))

	# limit of index for 1-hour running means
	indx_lim1 = sum(timehr <= (timehr[-1]-1.))

	# limit of index for 8-hour running means
	indx_lim8 = sum(timehr <= (timehr[-1]-8.))
	
	# prepare results matrices
	sim_PMc = np.zeros((indx_lim24))
	sim_PMf = np.zeros((indx_lim24))
	sim_SO2 = np.zeros((indx_lim15))
	sim_NO2 = np.zeros((indx_lim1))
	sim_O3 = np.zeros((indx_lim8))
	
	# prepare output concentrations (ug/m3)
	# gas-phase concentrations as ppb
	# index of O3, NO2, SO2
	O3_indx = comp_names.index('O3')
	NO2_indx = comp_names.index('NO2')
	SO2_indx = comp_names.index('SO2')
	
	# concentrations (# molecules/cm3)
	O3 = yrec[:, O3_indx]*Cfac
	NO2 = yrec[:, NO2_indx]*Cfac
	SO2 = yrec[:, SO2_indx]*Cfac
	
	# convert concentrations into ug/m3
	O3 = (O3/si.N_A)*1.e12*y_MW[0, O3_indx]
	NO2 = (NO2/si.N_A)*1.e12*y_MW[0, NO2_indx]
	SO2 = (SO2/si.N_A)*1.e12*y_MW[0, SO2_indx]
	
	PMf_indx = np.zeros((len(timehr), rbou_rec.shape[1]-1))
	PMc_indx = np.zeros((len(timehr), rbou_rec.shape[1]-1))

	# empty array for holding particulate matter mass concentration
	# PM2.5
	sim_PMf = np.zeros((len(timehr)))
	# PM10.0
	sim_PMc = np.zeros((len(timehr)))

	# get volumes in all particle size bins at all times (um3)
	vols = (4./3.*np.pi*cen_size**3.)
	
	# loop through times to get the indices for PM2.5 and PM10
	for it in range(len(timehr)):
		PMf_indxn = rbou_rec[it, 1::]*2. < 2.5
		PMc_indxn = rbou_rec[it, 1::]*2. < 10.

		# get the PM2.5 (ug/cm3) from simulated particle phase,
		# assuming a particle density of 1.4 g/cm3 (1.4e-6 ug/um3)
		sim_PMf[it] = (np.sum(Ndry[it, PMf_indxn]*vols[it, PMf_indxn]))*1.e-6

		# get the PM10. (ug/cm3) from simulated particle phase,
		# assuming a particle density of 1.4 g/cm3 (1.4e-6 ug/um3)
		sim_PMc[it] = (np.sum(Ndry[it, PMc_indxn]*vols[it, PMc_indxn]))*1.e-6
	
	# finally convert ug/cm3 to ug/m3
	sim_PMf = sim_PMf*1.e6
	sim_PMc = sim_PMc*1.e6

	# prepare for holding running means
	sim_PMf_mean = np.zeros((indx_lim24))
	sim_PMc_mean = np.zeros((indx_lim24))
	SO2_mean = np.zeros((indx_lim15))
	NO2_mean = np.zeros((indx_lim1))
	O3_mean = np.zeros((indx_lim8))

	# prepare for holding time relevant to running means
	time24 = np.zeros((indx_lim24))
	time15 = np.zeros((indx_lim15))
	time1 = np.zeros((indx_lim1))
	time8 = np.zeros((indx_lim8))

	# do running means
	for i in range(0, len(timehr[0:indx_lim15])):

		if (i < indx_lim24): # 24-hour running mean

			# get index of range to use
			tindx = (timehr>=timehr[i])*(timehr<=timehr[i]+24.)
			# convert to boolean
			tindx = (tindx==1)
			# estimate PM10 24-hour running mean (ug/m3)
			sim_PMc_mean[i] = sum(sim_PMc[tindx])/sum(tindx)
			# estimate PM2.5 24-hour running mean (ug/m3)
			sim_PMf_mean[i] = sum(sim_PMf[tindx])/sum(tindx)

			# get the corresponding time points (hour)
			time24[i] = np.mean(timehr[tindx])

		if (i < indx_lim15): # 15-minute running mean

			# get index of range to use
			tindx = (timemin>=timemin[i])*(timemin<=timemin[i]+15.)
			# convert to boolean
			tindx = (tindx==1)
			# estimate SO2 15-minute running mean (ug/m3)
			print(type(tindx))
			print(type(SO2))
			print(SO2_mean[i])
			print(sum(SO2[tindx]))
			print(sum(tindx))
			SO2_mean[i] = sum(SO2[tindx])/sum(tindx)
	
			# get the corresponding time points (hour)
			time15[i] = np.mean(timehr[tindx])
		
		if (i < indx_lim1): # 1-hour running mean

			# get index of range to use
			tindx = (timehr>=timehr[i])*(timehr<=timehr[i]+1.)
			# convert to boolean
			tindx = (tindx==1)
			# estimate NO2 1-hour running mean (ug/m3)
			NO2_mean[i] = sum(NO2[tindx])/sum(tindx)

			# get the corresponding time points (hour)
			time1[i] = np.mean(timehr[tindx])
		
		if (i < indx_lim8): # 8-hour running mean

			# get index of range to use
			tindx = (timehr>=timehr[i])*(timehr<=timehr[i]+8.)
			# convert to boolean
			tindx = (tindx==1)
			# estimate O3 8-hour running mean (ug/m3)
			O3_mean[i] = sum(O3[tindx])/sum(tindx)

			# get the corresponding time points (hour)
			time8[i] = np.mean(timehr[tindx])
	
	# PM10.0 24-hour running mean values (ug/m3)
	PMc_levels = [0., 16., 33., 50., 58., 66., 75., 83., 91., 100.]

	# PM2.5 24-hour running mean values
	PMf_levels = [0., 11., 23., 35., 41., 47., 53., 58., 64., 70.]

	# SO2 15-minute mean values
	SO2_levels = [0., 88., 177., 266., 354., 443., 532., 710., 887., 1064.]

	# NO2 hourly mean values
	NO2_levels = [0., 67., 134., 200., 267., 334., 400., 467., 534., 600.]

	# O3 8-hourly mean values
	O3_levels = [0., 33., 66., 100., 120., 140., 160., 187., 213., 240.]
	
	AQI_res = np.zeros((indx_lim24)) # prepare results matrix
	 
	# loop through times to get air quality index at each time
	for it in range(indx_lim24):

		# relevant time
		timen = time24[it]

		# find the index of 8 hour, 15 minunte and 1 hour means that
		# is closest to this 24 hour time
		indx8 = (np.abs(time8-timen) == (min(np.abs(time8-timen)))) == 1
		# in case more than one, choose the first
		indx8_first = np.where(indx8 == 1)[0][0]
		indx8[:] = 0; indx8[indx8_first] = True

		indx15 = (np.abs(time15-timen) == (min(np.abs(time15-timen)))) == 1
		# in case more than one, choose the first
		indx15_first = np.where(indx15 == 1)[0][0]
		indx15[:] = 0; indx15[indx15_first] = True

		indx1 = (np.abs(time1-timen) == (min(np.abs(time1-timen)))) == 1
		# in case more than one, choose the first
		indx1_first = np.where(indx1 == 1)[0][0]
		indx1[:] = 0; indx1[indx1_first] = True
		

		# get levels
		lev = []
		lev.append(np.sum(PMc_levels<=sim_PMc_mean[it]))
		lev.append(np.sum(PMf_levels<=sim_PMf_mean[it]))
		lev.append(np.sum(SO2_levels<=SO2_mean[indx15]))
		lev.append(np.sum(NO2_levels<=NO2_mean[indx1]))
		lev.append(np.sum(O3_levels<=O3_mean[indx8]))
		AQI_res[it] = max(lev)

	# plot air quality index as a function of time
	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	ax0.plot(time24[0: indx_lim24], AQI_res)

	ax0.set_ylabel(r'Air Quality Index', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	return() # end function

	# define function to plot total volatile organic 
	# components in gas-phase against time
def tvoc_calc(self):

	# inputs: ------------------------------------------------------------------
	# self - reference to GUI
	# --------------------------------------------------------------------------
	
	# get required variables from self
	yrec = np.zeros((self.ro_obj.yrec.shape[0], self.ro_obj.yrec.shape[1]))
	yrec[:, :] = self.ro_obj.yrec[:, :]
	num_comp = self.ro_obj.nc
	num_sb = self.ro_obj.nsb
	timehr = self.ro_obj.thr
	comp_names = self.ro_obj.names_of_comp
	rel_SMILES = self.ro_obj.rSMILES
	y_MW = (np.array((self.ro_obj.comp_MW))).reshape(1, -1)
	H2Oi = self.ro_obj.H2O_ind
	seedi = self.ro_obj.seed_ind
	indx_plot = self.ro_obj.plot_indx
	comp0 = self.ro_obj.init_comp
	rbou_rec = np.zeros((self.ro_obj.rad.shape[0], self.ro_obj.rad.shape[1]))
	rbou_rec[:, :] = self.ro_obj.rad[:, :]
	space_mode = self.ro_obj.spacing
	HC = np.array((self.ro_obj.HyC))

	# get index of components that are hydrocarbons
	HCindx = np.array(((HC > 0.)))
	# prepare for indexing y
	HCindx = np.squeeze(HCindx.reshape(1, -1))

	HCsum = np.zeros((len(timehr), 1))

	for ti in range(len(timehr)):
		# get relevant concentrations in gas-phase (ppb)
		HCsum[ti, 0] = np.sum(yrec[ti, 0:num_comp][HCindx])

	# plot against time
	plt.ion() # show results to screen and turn on interactive mode
		
	# prepare plot
	fig, (ax0) = plt.subplots(1, 1, figsize=(14, 7))
	ax0.plot(timehr, HCsum)

	ax0.set_ylabel(r'TVOC (inc. CH4) (ppb)', fontsize = 14)
	ax0.set_xlabel(r'Time through simulation (hours)', fontsize = 14)
	ax0.yaxis.set_tick_params(labelsize = 14, direction = 'in')
	ax0.xaxis.set_tick_params(labelsize = 14, direction = 'in')

	return() # end function

# function to makes spines invisible
def make_patch_spines_invisible(ax):
	ax.set_frame_on(True)
	ax.patch.set_visible(False)
	for sp in ax.spines.values():
		sp.set_visible(False)

# function for doing colorbar tick labels in standard notation
def fmt(x, pos):
	a, b = '{:.1e}'.format(x).split('e')
	b = int(b)
	return r'${} \times 10^{{{}}}$'.format(a, b)
