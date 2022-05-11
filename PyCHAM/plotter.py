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

	# chamber condition ---------------------------------------------------------
	# retrieve results
	(num_sb, num_comp, Cfac, yrec, Ndry, rbou_rec, x, timehr, _, 
		y_mw, Nwet, _, y_MV, _, wall_on, space_mode, indx_plot, 
		comp0, _, PsatPa, OC, _, _, _, _, _, _, _) = retr_out.retr_out(dir_path)
	
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
		if not (indx_plot):
			mess = str('Please note, no initial gas-phase concentrations were registered, therefore the gas-phase standard plot will not be shown')
			self.l203a.setText(mess)
			# if there are no gaseous components then prepare figure
			fig, (ax1) = plt.subplots(1, 1, figsize=(14, 7))
		else:
			# if there are both gaseous components and particle size bins then prepare figure
			fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(14, 7))

		# parasite axis setup --------------------------------------------------------------
		par1 = ax1.twinx() # first parasite axis
		par2 = ax1.twinx() # second parasite axis
		
		# Offset the right spine of par2.  The ticks and label have already been
		# placed on the right by twinx above.
		par2.spines["right"].set_position(("axes", 1.2))
		# Having been created by twinx, par2 has its frame off, so the line of its
		# detached spine is invisible.  First, activate the frame but make the patch
		# and spines invisible.
		make_patch_spines_invisible(par2)
		# second, show the right spine
		par2.spines["right"].set_visible(True)	
		# ----------------------------------------------------------------------------------------

	if (indx_plot):
		
		
		ymax = 0. # start tracking maximum value for plot label
		
		# action units for gas-phase concentrations
		if (uc == 0): # ppb
			gp_conc = yrec[:, 0:num_comp] # ppb is original units
			gpunit = '(ppb)'
		if (uc == 1 or uc == 2): # ug/m3 or # molecules/cm3

			y_MW = np.array(y_mw) # convert to numpy array from list
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
			ax0.text(x = timehr[0]-(timehr[-1]-timehr[0])/9.5, y = ymax*1.05, s='a)', size=14)

		# end of gas-phase concentration sub-plot ---------------------------------------
	
	# particle properties sub-plot --------------------------------------------------
	if (num_asb > 0): # if particles present		

		if (timehr.ndim == 0): # occurs if only one time step saved
			Ndry = np.array(Ndry.reshape(1, num_asb))
			x = np.array(x.reshape(1, num_asb))
			rbou_rec = np.array(rbou_rec.reshape(1, num_sb))
		if (num_asb == 1): # just one particle size bin (wall included in num_sb)
			Ndry = np.array(Ndry.reshape(len(timehr), num_asb))
			x = np.array(x.reshape(len(timehr), num_asb))

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
	
		# set contour levels
		levels = (MaxNLocator(nbins = 100).tick_values(np.min(z[~np.isnan(z)]), 
				np.max(z[~np.isnan(z)])))
	
		# associate colours and contour levels
		norm1 = BoundaryNorm(levels, ncolors=cm.N, clip=True)
		
		# contour plot with times (hours) along x axis and 
		# particle diameters (nm) along y axis
		for ti in range(len(timehr)-1): # loop through times
			p1 = ax1.pcolormesh(timehr[ti:ti+2], (rbou_rec[ti, :]*2*1e3), z[:, ti].reshape(-1, 1), cmap=cm, norm=norm1)
	
		# if logarithmic spacing of size bins specified, plot vertical axis 
		# logarithmically
		if space_mode == 'log':
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
		
		cb = plt.colorbar(p1, format=ticker.FuncFormatter(fmt), pad=0.25, ax=ax1)
		cb.ax.tick_params(labelsize=14)   
		# colour bar label
		cb.set_label('dN (#$\,$$\mathrm{cm^{-3}}$)/d$\,$log$_{10}$(D$\mathrm{_p}$ ($\mathrm{\mu m}$))', size=14, rotation=270, labelpad=20)

		# ----------------------------------------------------------------------------------------
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

		# mass concentration of particles ---------------------------------------------------------------
		# array for mass concentration with time
		MCvst = np.zeros((1, len(timehr)))
		
		# first obtain just the particle-phase concentrations (molecules/cm3)
		yrp = yrec[:, num_comp:num_comp*(num_asb+1)]
		# loop through size bins to convert to ug/m3
		for sbi in range(num_asb):
			yrp[:, sbi*num_comp:(sbi+1)*num_comp] = ((yrp[:, sbi*num_comp:(sbi+1)*num_comp] /si.N_A)*y_mw)*1.e12
		
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
